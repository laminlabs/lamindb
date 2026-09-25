"""Sync Notion pages to LaminDB records.

.. autofunction:: sync_from_notion
.. autoclass:: SyncReport

"""

from __future__ import annotations

import json
import os
import re
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any
from urllib.parse import quote, urlparse

import httpx
from lamin_utils import logger
from rich.console import Console
from rich.markup import escape as rich_escape

import lamindb as ln

from ..base.types import PROJECT_STATUS_TO_CODE as LAMIN_PROJECT_STATUS_TO_CODE

API_VERSION = "2026-03-11"
BASE = "https://api.notion.com/v1"
UUID_DASHED_PATTERN = re.compile(
    r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
)
UUID_COMPACT_PATTERN = re.compile(r"^[0-9a-fA-F]{32}$")
MARKDOWN_IMAGE_LINK_PATTERN = re.compile(r"!\[([^\]]*)\]\((https?://[^)\s]+)\)")
HTML_IMAGE_SRC_PATTERN = re.compile(
    r'(<img\b[^>]*\bsrc=["\'])(https?://[^"\']+)(["\'][^>]*>)',
    flags=re.IGNORECASE,
)
RICH_CONSOLE = Console(force_terminal=True, no_color=False)
NOTION_EMBEDDED_IMAGE_WIDTH_PX = 500


def _compact_uuid(value: str) -> str:
    """Format UUID-like strings without dashes for CLI-facing messages."""
    return value.replace("-", "") if UUID_DASHED_PATTERN.match(value) else value


def _normalize_notion_id(value: str | None) -> str | None:
    if not isinstance(value, str):
        return None
    value = value.strip()
    if not value:
        return None
    return _compact_uuid(value)


def _notion_api_id(value: str) -> str:
    """Canonical dashed UUID form for Notion API path parameters."""
    if UUID_COMPACT_PATTERN.match(value):
        return (
            f"{value[0:8]}-{value[8:12]}-{value[12:16]}-{value[16:20]}-{value[20:32]}"
        )
    return value


def _extract_emoji(payload: dict) -> str | None:
    icon = payload.get("icon")
    if not isinstance(icon, dict) or icon.get("type") != "emoji":
        return None
    emoji = icon.get("emoji")
    if not isinstance(emoji, str):
        return None
    emoji = emoji.strip()
    return emoji or None


@dataclass
class SyncReport:
    apply: bool = False
    message: str | None = None
    discovered_pages: int = 0
    discovered: int = 0
    created: int = 0
    updated: int = 0
    unchanged: int = 0
    pending_relations: int = 0
    failed: int = 0
    errors: list[str] = field(default_factory=list)
    databases: list[str] = field(default_factory=list)
    created_record_types: list[str] = field(default_factory=list)
    create_record_types: list[str] = field(default_factory=list)
    updated_record_types: list[str] = field(default_factory=list)
    update_record_types: list[str] = field(default_factory=list)
    created_feature_types: list[str] = field(default_factory=list)
    create_feature_types: list[str] = field(default_factory=list)
    created_ulabel_types: list[str] = field(default_factory=list)
    create_ulabel_types: list[str] = field(default_factory=list)
    created_ulabels: list[str] = field(default_factory=list)
    create_ulabels: list[str] = field(default_factory=list)
    created_schemas: list[str] = field(default_factory=list)
    create_schemas: list[str] = field(default_factory=list)
    updated_schemas: list[str] = field(default_factory=list)
    update_schemas: list[str] = field(default_factory=list)
    created_features: list[str] = field(default_factory=list)
    create_features: list[str] = field(default_factory=list)
    updated_features: list[str] = field(default_factory=list)
    update_features: list[str] = field(default_factory=list)
    created_artifacts: list[str] = field(default_factory=list)
    create_artifacts: list[str] = field(default_factory=list)
    created_relation_stubs: list[str] = field(default_factory=list)
    create_relation_stubs: list[str] = field(default_factory=list)
    relation_value_links: list[str] = field(default_factory=list)
    mapped_project_record_relations: list[str] = field(default_factory=list)
    created_projects: int = 0
    updated_projects: int = 0
    unchanged_projects: int = 0
    unmapped_properties: list[str] = field(default_factory=list)

    def to_pretty_text(self) -> str:
        """Render a concise human-readable sync report."""

        def safe(value: str | int) -> str:
            return rich_escape(str(value))

        def metric(key: str, value: str | int, color: str = "white") -> str:
            return f"[bold white]{key}[/]: [{color}]{safe(value)}[/]"

        lines: list[str] = []
        if not self.apply:
            lines.append(
                "[bold yellow]Dry run: nothing got created. If you're happy, pass apply=True or --apply on the CLI.[/]"
            )
        else:
            lines.append("[bold cyan]Sync report[/]")
        lines.append("")
        lines.append("[bold cyan]Scope[/]")
        if self.discovered_pages > 0:
            lines.append(f"[bold]Discovered {self.discovered_pages} Notion pages.[/]")

        lines.append("")
        lines.append("[bold cyan]Actions[/]")
        action_color = "yellow" if not self.apply else "green"
        if self.create_record_types:
            lines.append(
                metric(
                    "create_record_types",
                    ", ".join(self.create_record_types),
                    action_color,
                )
            )
        if self.created_record_types:
            lines.append(
                metric(
                    "created_record_types",
                    ", ".join(self.created_record_types),
                    "green",
                )
            )
        if self.update_record_types:
            lines.append(
                metric(
                    "update_record_types",
                    ", ".join(self.update_record_types),
                    action_color,
                )
            )
        if self.updated_record_types:
            lines.append(
                metric(
                    "updated_record_types",
                    ", ".join(self.updated_record_types),
                    "green",
                )
            )
        if self.create_feature_types:
            lines.append(
                metric(
                    "create_feature_types",
                    ", ".join(self.create_feature_types),
                    action_color,
                )
            )
        if self.created_feature_types:
            lines.append(
                metric(
                    "created_feature_types",
                    ", ".join(self.created_feature_types),
                    "green",
                )
            )
        if self.create_ulabel_types:
            lines.append(
                metric(
                    "create_ulabel_types",
                    ", ".join(self.create_ulabel_types),
                    action_color,
                )
            )
        if self.created_ulabel_types:
            lines.append(
                metric(
                    "created_ulabel_types",
                    ", ".join(self.created_ulabel_types),
                    "green",
                )
            )
        if self.create_schemas:
            lines.append(
                metric("create_schemas", ", ".join(self.create_schemas), action_color)
            )
        if self.created_schemas:
            lines.append(
                metric("created_schemas", ", ".join(self.created_schemas), "green")
            )
        if self.update_schemas:
            lines.append(
                metric("update_schemas", ", ".join(self.update_schemas), action_color)
            )
        if self.updated_schemas:
            lines.append(
                metric("updated_schemas", ", ".join(self.updated_schemas), "green")
            )
        if self.create_features:
            lines.append("[bold]create_features[/]:")
            lines.extend(
                f"  [{action_color}]{safe(feature)}[/]"
                for feature in self.create_features
            )
        if self.created_features:
            lines.append("[bold]created_features[/]:")
            lines.extend(
                f"  [green]{safe(feature)}[/]" for feature in self.created_features
            )
        if self.update_features:
            lines.append("[bold]update_features[/]:")
            lines.extend(
                f"  [{action_color}]{safe(feature)}[/]"
                for feature in self.update_features
            )
        if self.updated_features:
            lines.append("[bold]updated_features[/]:")
            lines.extend(
                f"  [green]{safe(feature)}[/]" for feature in self.updated_features
            )
        if self.create_ulabels:
            lines.append("[bold]create_ulabels[/]:")
            lines.extend(
                f"  [{action_color}]{safe(label)}[/]" for label in self.create_ulabels
            )
        if self.created_ulabels:
            lines.append("[bold]created_ulabels[/]:")
            lines.extend(f"  [green]{safe(label)}[/]" for label in self.created_ulabels)
        if self.create_artifacts:
            lines.append("[bold]create_artifacts[/]:")
            lines.extend(
                f"  [{action_color}]{safe(artifact)}[/]"
                for artifact in self.create_artifacts
            )
        if self.created_artifacts:
            lines.append("[bold]created_artifacts[/]:")
            lines.extend(
                f"  [green]{safe(artifact)}[/]" for artifact in self.created_artifacts
            )
        if self.create_relation_stubs:
            lines.append("[bold]create_relation_stubs[/]:")
            lines.extend(
                f"  [{action_color}]{safe(stub)}[/]"
                for stub in self.create_relation_stubs
            )
        if self.created_relation_stubs:
            lines.append("[bold]created_relation_stubs[/]:")
            lines.extend(
                f"  [green]{safe(stub)}[/]" for stub in self.created_relation_stubs
            )
        if self.relation_value_links:
            lines.append("[bold]relation_value_links[/]:")
            lines.extend(
                f"  [{action_color}]{safe(summary)}[/]"
                for summary in self.relation_value_links
            )
        if self.mapped_project_record_relations:
            lines.append("[bold]mapped_project_record_relations[/]:")
            lines.extend(
                f"  [{action_color}]{safe(detail)}[/]"
                for detail in self.mapped_project_record_relations
            )
        if self.unmapped_properties:
            lines.append("[bold]unmapped_properties[/]:")
            lines.extend(
                f"  [yellow]{safe(detail)}[/]" for detail in self.unmapped_properties
            )
        lines.extend(
            [
                metric("create_records", self.created, action_color),
                metric("update_records", self.updated, action_color),
                metric("unchanged_records", self.unchanged, action_color),
                metric("create_projects", self.created_projects, action_color),
                metric("update_projects", self.updated_projects, action_color),
                metric("unchanged_projects", self.unchanged_projects, action_color),
                metric("pending_relations", self.pending_relations, action_color),
                metric("failed_records", self.failed, action_color),
            ]
        )
        if self.errors:
            lines.append("")
            lines.append("[bold red]Errors[/]")
            lines.extend(f"[red]- {error}[/]" for error in self.errors)
        return "\n".join(lines)


def _flatten(prop: dict, *, people_name_cache: dict[str, str] | None = None) -> Any:
    t = prop.get("type")
    if t in ("title", "rich_text"):
        return "".join(s.get("plain_text", "") for s in (prop.get(t) or []))
    if t in (
        "email",
        "phone_number",
        "url",
        "number",
        "checkbox",
        "created_time",
        "last_edited_time",
    ):
        return prop.get(t)
    if t in ("select", "status"):
        opt = prop.get(t)
        return opt["name"] if opt else None
    if t == "multi_select":
        return [o["name"] for o in (prop.get(t) or [])]
    if t == "date":
        d = prop.get("date")
        return d["start"] if d else None
    if t in ("people", "relation"):
        ids: list[str] = []
        for item in prop.get(t) or []:
            raw_id = item.get("id")
            if isinstance(raw_id, str) and raw_id:
                normalized_id = _normalize_notion_id(raw_id) or raw_id
                ids.append(normalized_id)
                if t == "people" and people_name_cache is not None:
                    raw_name = item.get("name")
                    if isinstance(raw_name, str):
                        person_name = raw_name.strip()
                        if person_name:
                            people_name_cache[normalized_id] = person_name
        return ids
    if t in ("created_by", "last_edited_by"):
        u = prop.get(t)
        user_id = u.get("id") if isinstance(u, dict) else None
        if not isinstance(user_id, str) or not user_id:
            return None
        return _normalize_notion_id(user_id) or user_id
    if t == "files":
        out = []
        for f in prop.get("files") or []:
            src = f.get(f.get("type"), {})
            if "url" in src:
                out.append(src["url"])
        return out
    if t == "formula":
        formula = prop.get("formula")
        if not isinstance(formula, dict):
            return None
        formula_type = formula.get("type")
        if formula_type in {"string", "number", "boolean"}:
            return formula.get(formula_type)
        if formula_type == "date":
            date_value = formula.get("date")
            if isinstance(date_value, dict):
                return date_value.get("start")
            return None
        return None
    return None  # rollup, formula, unknown


def _page_title(page: dict) -> str:
    """The title of a page, whatever the title property happens to be called."""
    for prop in page.get("properties", {}).values():
        if prop.get("type") == "title":
            return _flatten(prop)
    return ""


def _parse_notion_timestamp(value: Any) -> datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None


def _normalized_timestamp(value: datetime | None) -> datetime | None:
    if value is None:
        return None
    if value.tzinfo is None:
        value = value.replace(tzinfo=UTC)
    return value.astimezone(UTC).replace(microsecond=0)


class _NotionReader:
    """Read-only Notion reader. Databases contain data sources; rows live on the data source."""

    def __init__(self, token: str) -> None:
        if not token:
            raise ValueError("A Notion access token is required.")
        self.s = httpx.Client(
            transport=httpx.HTTPTransport(verify=True, http2=False, trust_env=True)
        )
        self.s.headers.update(
            {
                "Authorization": f"Bearer {token}",
                "Notion-Version": API_VERSION,
                "Content-Type": "application/json",
            }
        )
        self._ds: dict[str, str] = {}
        self._schema: dict[str, dict] = {}
        self._titles: dict[str, dict[str, str]] = {}
        self._notion_people_names_by_id: dict[str, str] = {}

    def _call(
        self,
        method: str,
        path: str,
        body: dict | None = None,
        params: dict | None = None,
    ) -> dict:
        # Notion caps at ~3 req/s and returns 429 with Retry-After. Honour it,
        # and retry transient 5xx with backoff, so a real workspace doesn't die.
        for attempt in range(6):
            r = self.s.request(
                method, f"{BASE}{path}", json=body, params=params, timeout=30
            )
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", "1")))
                continue
            if r.status_code == 401:
                raise PermissionError("Invalid or expired Notion token.")
            if r.status_code == 404:
                raise LookupError(
                    f"404 on {path} — not found, or not shared with this connection "
                    "(Notion: ••• -> Connections -> Connect to)."
                )
            if r.status_code >= 500:
                time.sleep(2**attempt)
                continue
            r.raise_for_status()
            return r.json()
        raise RuntimeError(f"Notion API: gave up after retries on {path}")

    def data_sources(self, database_id: str) -> list[dict]:
        """List the data sources under a database: [{id, name}, ...]."""
        return self._call("GET", f"/databases/{database_id}").get("data_sources", [])

    def _resolve(self, database_id: str) -> str:
        if database_id not in self._ds:
            sources = self.data_sources(database_id)
            if not sources:
                raise LookupError(
                    f"No data sources on database {_compact_uuid(database_id)!r}."
                )
            if len(sources) > 1:
                logger.warning(
                    f"database {_compact_uuid(database_id)!r} has {len(sources)} data sources; "
                    f"using {sources[0]['name']!r}"
                )
            self._ds[database_id] = sources[0]["id"]
        return self._ds[database_id]

    def schema(self, database_id: str) -> dict[str, dict]:
        """{property_name: {"type": str, "target": str | None, "dual": dict | None, "choices": list[str] | None, "formula_expression": str | None}}.

        `target` and `dual` are set only for relation properties. `target` names
        the data source the relation points at. `dual` is Notion's synced-property
        info when the relation is two-way — the two sides describe the same edge,
        so a sync should follow only one of them. Cached.
        """
        if database_id in self._schema:
            return self._schema[database_id]
        ds = self._resolve(database_id)
        props = self._call("GET", f"/data_sources/{ds}").get("properties", {})
        out: dict[str, dict] = {}
        for name, p in props.items():
            t = p.get("type", "")
            target = None
            dual = None
            choices = None
            formula_expression = None
            if t == "formula":
                formula = p.get("formula")
                if isinstance(formula, dict):
                    expr = formula.get("expression")
                    if isinstance(expr, str) and expr.strip():
                        formula_expression = expr.strip()
            if t == "relation":
                rel = p.get("relation", {})
                target = rel.get("data_source_id") or rel.get("database_id")
                dual = rel.get("dual_property")
            if t in {"select", "status", "multi_select"}:
                type_payload = p.get(t, {})
                options = (
                    type_payload.get("options", [])
                    if isinstance(type_payload, dict)
                    else []
                )
                choices = [
                    option.get("name", "").strip()
                    for option in options
                    if isinstance(option, dict) and isinstance(option.get("name"), str)
                ]
                choices = [value for value in choices if value]
            out[name] = {
                "type": t,
                "target": target,
                "dual": dual,
                "choices": choices,
                "formula_expression": formula_expression,
            }
        self._schema[database_id] = out
        return out

    def columns(self, database_id: str) -> dict[str, str]:
        """Return {property_name: notion_type}."""
        return {k: v["type"] for k, v in self.schema(database_id).items()}

    def _query(self, ds: str, limit: int | None = None) -> list[dict]:
        """Pages in a data source, raw. Paginates until exhausted or `limit` reached."""
        if limit is not None and limit <= 0:
            return []
        page_size = 100 if limit is None else min(100, max(1, limit))
        body: dict[str, Any] = {"page_size": page_size}
        pages: list[dict] = []
        while True:
            payload = self._call("POST", f"/data_sources/{ds}/query", body)
            pages.extend(payload.get("results", []))
            if limit is not None and len(pages) >= limit:
                return pages[:limit]
            if not payload.get("has_more"):
                return pages
            body["start_cursor"] = payload["next_cursor"]

    def title_map(self, ds_id: str) -> dict[str, str]:
        """{page_id: title} for a data source. One query per source, cached."""
        if ds_id not in self._titles:
            try:
                self._titles[ds_id] = {
                    p["id"]: _page_title(p) for p in self._query(ds_id)
                }
            except LookupError:
                logger.warning(
                    f"data source {ds_id} is not shared with this connection; "
                    "relation IDs pointing at it cannot be titled"
                )
                self._titles[ds_id] = {}
        return self._titles[ds_id]

    def rows(
        self,
        database_id: str,
        drop: set[str] | None = None,
        limit: int | None = None,
        *,
        include_page_emoji: bool = False,
    ) -> list[dict]:
        """Every page flattened to a dict.

        Relation values are always lists of Notion page UUIDs — never titles.
        UUIDs are the stable join key; titles are not unique and change on
        rename. Use :meth:`relation_titles` for display names.

        Args:
            database_id: the Notion database ID.
            drop: property names or property types to omit, e.g.
                {"messages", "created_by"}.
            limit: stop after this many rows. Useful for sampling a large
                database without paginating through all of it.
        """
        drop = drop or set()
        ds = self._resolve(database_id)

        rows: list[dict] = []
        for page in self._query(ds, limit=limit):
            row: dict[str, Any] = {
                "notion_id": None,
                "created_time": None,
                "last_edited_time": None,
            }
            for name, prop in page.get("properties", {}).items():
                if name in drop or prop.get("type") in drop:
                    continue
                row[name] = _flatten(
                    prop, people_name_cache=self._notion_people_names_by_id
                )
            # page-level fields win over any same-named user property
            page_id = page.get("id")
            row["notion_id"] = (
                _normalize_notion_id(page_id) if isinstance(page_id, str) else page_id
            )
            row["created_time"] = page.get("created_time")
            row["last_edited_time"] = page.get("last_edited_time")
            if include_page_emoji:
                row["__notion_emoji__"] = _extract_emoji(page)
            rows.append(row)
        return rows

    @staticmethod
    def _rich_text_to_markdown(rich_text: list[dict]) -> str:
        chunks: list[str] = []
        for segment in rich_text:
            text = segment.get("plain_text", "")
            if not text:
                continue
            href = segment.get("href")
            if isinstance(href, str) and href and href != text:
                text = f"[{text}]({href})"
            annotations = segment.get("annotations") or {}
            if annotations.get("code"):
                text = f"`{text}`"
            if annotations.get("bold"):
                text = f"**{text}**"
            if annotations.get("italic"):
                text = f"*{text}*"
            if annotations.get("strikethrough"):
                text = f"~~{text}~~"
            chunks.append(text)
        return "".join(chunks)

    @staticmethod
    def _block_file_url(payload: dict | None) -> str | None:
        if not isinstance(payload, dict):
            return None
        source_type = payload.get("type")
        if not isinstance(source_type, str):
            return None
        source_payload = payload.get(source_type)
        if not isinstance(source_payload, dict):
            return None
        url = source_payload.get("url")
        if not isinstance(url, str):
            return None
        url = url.strip()
        return url or None

    def _block_caption_markdown(self, payload: dict | None) -> str:
        if not isinstance(payload, dict):
            return ""
        caption = payload.get("caption") or []
        if not isinstance(caption, list):
            return ""
        return self._rich_text_to_markdown(caption).strip()

    @staticmethod
    def _table_cell_markdown(cell: Any) -> str:
        if not isinstance(cell, list):
            return ""
        text = _NotionReader._rich_text_to_markdown(cell).strip()
        text = text.replace("|", r"\|")
        # Preserve line breaks in markdown table cells.
        return "<br>".join(text.splitlines()) if text else ""

    def _table_to_markdown_lines(
        self, table_payload: dict | None, table_rows: list[dict], indent: str
    ) -> list[str]:
        if not table_rows:
            return []

        parsed_rows: list[list[str]] = []
        max_cols = 0
        for row in table_rows:
            row_payload = row.get("table_row", {}) if isinstance(row, dict) else {}
            cells = (
                row_payload.get("cells", []) if isinstance(row_payload, dict) else []
            )
            rendered = [self._table_cell_markdown(cell) for cell in cells]
            parsed_rows.append(rendered)
            max_cols = max(max_cols, len(rendered))

        if max_cols == 0:
            return []

        normalized_rows = [row + [""] * (max_cols - len(row)) for row in parsed_rows]
        has_header = bool(
            isinstance(table_payload, dict) and table_payload.get("has_column_header")
        )
        if has_header:
            header = normalized_rows[0]
            body = normalized_rows[1:]
        else:
            header = normalized_rows[0]
            body = normalized_rows[1:]

        def render_row(values: list[str]) -> str:
            return f"{indent}| " + " | ".join(values) + " |"

        lines = [render_row(header), render_row(["---"] * max_cols)]
        lines.extend(render_row(values) for values in body)
        return lines

    def _block_to_markdown_lines(
        self, block: dict, *, depth: int = 0, parent_is_numbered: bool = False
    ) -> list[str]:
        block_type = block.get("type")
        payload = block.get(block_type) if isinstance(block_type, str) else None
        rich_text = payload.get("rich_text", []) if isinstance(payload, dict) else []
        text = self._rich_text_to_markdown(rich_text)
        indent = "  " * depth

        lines: list[str] = []
        if block_type == "heading_1":
            lines.append(f"{indent}# {text}".rstrip())
        elif block_type == "heading_2":
            lines.append(f"{indent}## {text}".rstrip())
        elif block_type == "heading_3":
            lines.append(f"{indent}### {text}".rstrip())
        elif block_type == "paragraph":
            if text:
                lines.append(f"{indent}{text}".rstrip())
                # Notion paragraph blocks map to separate markdown paragraphs.
                # Keep a blank line separator so adjacent paragraph blocks don't collapse.
                lines.append("")
        elif block_type == "bulleted_list_item":
            lines.append(f"{indent}- {text}".rstrip())
        elif block_type == "numbered_list_item":
            lines.append(f"{indent}1. {text}".rstrip())
        elif block_type == "to_do":
            checked = (
                bool(payload.get("checked")) if isinstance(payload, dict) else False
            )
            marker = "x" if checked else " "
            lines.append(f"{indent}- [{marker}] {text}".rstrip())
        elif block_type == "quote":
            # Notion quote rich_text can contain embedded newlines; markdown
            # requires each line to be prefixed with ">" to stay in the quote.
            quote_lines = text.splitlines() or [""]
            for quote_line in quote_lines:
                lines.append(f"{indent}> {quote_line}".rstrip())
        elif block_type == "code":
            language = payload.get("language") if isinstance(payload, dict) else None
            fence = f"```{language}" if language else "```"
            lines.extend([f"{indent}{fence}", f"{indent}{text}", f"{indent}```"])
        elif block_type == "image":
            image_url = self._block_file_url(payload)
            if image_url is not None:
                caption = self._block_caption_markdown(payload) or "image"
                lines.append(f"{indent}![{caption}]({image_url})")
        elif block_type in {"file", "pdf", "video", "audio"}:
            file_url = self._block_file_url(payload)
            if file_url is not None:
                caption = self._block_caption_markdown(payload)
                label = caption or _short_file_source(file_url)
                lines.append(f"{indent}[{label}]({file_url})")
        elif block_type in {"bookmark", "embed", "link_preview"}:
            # Keep external previews (e.g. Slack URLs) as plain links in markdown.
            preview_url = payload.get("url") if isinstance(payload, dict) else None
            if isinstance(preview_url, str):
                preview_url = preview_url.strip()
                if preview_url:
                    lines.append(f"{indent}{preview_url}")
        elif block_type == "table_of_contents":
            lines.append(f"{indent}<!-- display-table-of-contents -->")
        elif block_type == "divider":
            lines.append(f"{indent}---")
        else:
            if text:
                prefix = f"{indent}1. " if parent_is_numbered else indent
                lines.append(f"{prefix}{text}".rstrip())

        if block_type == "toggle":
            lines.append(f"{indent}<details>")
            summary = text.strip() or "Details"
            lines.append(f"{indent}<summary>{summary}</summary>")
            lines.append(f"{indent}<p>")
            if block.get("has_children") and block.get("id"):
                child_blocks = self._iter_block_children(block["id"])
                for child in child_blocks:
                    lines.extend(
                        self._block_to_markdown_lines(
                            child,
                            depth=depth + 1,
                            parent_is_numbered=False,
                        )
                    )
            lines.append(f"{indent}</p>")
            lines.append(f"{indent}</details>")
            return lines

        if block_type == "table":
            table_rows = (
                self._iter_block_children(block["id"])
                if block.get("has_children") and block.get("id")
                else []
            )
            lines.extend(self._table_to_markdown_lines(payload, table_rows, indent))
            return lines

        if block.get("has_children") and block.get("id"):
            child_blocks = self._iter_block_children(block["id"])
            for child in child_blocks:
                lines.extend(
                    self._block_to_markdown_lines(
                        child,
                        depth=depth + 1,
                        parent_is_numbered=block_type == "numbered_list_item",
                    )
                )
        return lines

    def _iter_block_children(self, block_id: str) -> list[dict]:
        children: list[dict] = []
        cursor: str | None = None
        while True:
            params: dict[str, Any] = {"page_size": 100}
            if cursor:
                params["start_cursor"] = cursor
            payload = self._call(
                "GET", f"/blocks/{_notion_api_id(block_id)}/children", params=params
            )
            children.extend(payload.get("results", []))
            if not payload.get("has_more"):
                return children
            cursor = payload.get("next_cursor")

    def page_markdown(self, page_id: str) -> str:
        """Page body as markdown, including nested child blocks."""
        markdown_lines: list[str] = []
        for block in self._iter_block_children(page_id):
            markdown_lines.extend(self._block_to_markdown_lines(block))
        text = "\n".join(line.rstrip() for line in markdown_lines).strip()
        return re.sub(r"\n{3,}", "\n\n", text)

    def relation_titles(
        self,
        database_id: str,
        drop: set[str] | None = None,
    ) -> dict[str, str]:
        """{page_uuid: title} for every page reachable via this database's relations.

        A display-only side table: join it against the UUIDs in :meth:`rows`.
        Targets that are not shared with the connection contribute nothing and
        emit a warning.
        """
        drop = drop or set()
        out: dict[str, str] = {}
        for name, spec in self.schema(database_id).items():
            if spec["type"] != "relation" or not spec["target"] or name in drop:
                continue
            out.update(self.title_map(spec["target"]))
        return out

    def to_json(
        self,
        database_id: str,
        path: str | None = None,
        titles: bool = True,
        drop: set[str] | None = None,
    ) -> str:
        """Columns, rows, and a relation-title lookup as JSON.

        Rows carry relation UUIDs. `titles` maps those UUIDs to display names;
        set `titles=False` to skip the extra queries.
        """
        doc: dict[str, Any] = {
            "database_id": database_id,
            "columns": self.columns(database_id),
            "titles": self.relation_titles(database_id, drop=drop) if titles else {},
            "rows": self.rows(database_id, drop=drop),
        }
        text = json.dumps(doc, indent=2, ensure_ascii=False)
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
        return text

    def page_text(self, page_id: str) -> str:
        """Page body (meeting notes) as plain text. Top-level blocks only.

        Content is separate from properties; this reads block children and joins
        their rich_text. Nested/toggle children are not descended into.
        """
        lines: list[str] = []
        cursor: str | None = None
        while True:
            params: dict[str, Any] = {"page_size": 100}
            if cursor:
                params["start_cursor"] = cursor
            payload = self._call(
                "GET", f"/blocks/{_notion_api_id(page_id)}/children", params=params
            )
            for block in payload.get("results", []):
                t = block.get("type")
                data = block.get(t)
                rich = data.get("rich_text") if isinstance(data, dict) else None
                if rich:
                    text = "".join(s.get("plain_text", "") for s in rich)
                    if text:
                        lines.append(text)
            if not payload.get("has_more"):
                break
            cursor = payload.get("next_cursor")
        return "\n".join(lines)

    def page(self, page_id: str, drop: set[str] | None = None) -> dict:
        """Fetch one page by ID and flatten it to a dict, like a :meth:`rows` entry.

        Args:
            page_id: Notion page UUID.
            drop: property names or property types to omit, same semantics as
                :meth:`rows`.

        Returns:
            Dict with one key per property plus ``notion_id`` and
            page-level timestamps (``created_time``, ``last_edited_time``).
            Relation/people values are lists of Notion page UUIDs — the stable
            join key.
        """
        drop = drop or set()
        raw = self._call("GET", f"/pages/{_notion_api_id(page_id)}")
        row: dict[str, Any] = {
            "notion_id": None,
            "created_time": None,
            "last_edited_time": None,
        }
        for name, prop in raw.get("properties", {}).items():
            if name in drop or prop.get("type") in drop:
                continue
            row[name] = _flatten(
                prop, people_name_cache=self._notion_people_names_by_id
            )
        raw_id = raw.get("id")
        row["notion_id"] = (
            _normalize_notion_id(raw_id) if isinstance(raw_id, str) else raw_id
        )
        row["created_time"] = raw.get("created_time")
        row["last_edited_time"] = raw.get("last_edited_time")
        return row


@contextmanager
def _bulk_creation():
    """Disable the per-create similar-name search for the duration of a bulk sync."""
    prev = ln.settings.creation.search_names
    ln.settings.creation.search_names = False
    try:
        yield
    finally:
        ln.settings.creation.search_names = prev


def _kinds(spec: dict) -> tuple[set, set, set]:
    """Split a Notion schema into (relation-props, label-props, file-props)."""
    rel = {p for p, s in spec.items() if s["type"] in ("relation", "people")}
    lab = {
        p for p, s in spec.items() if s["type"] in ("select", "status", "multi_select")
    }
    file = {p for p, s in spec.items() if s["type"] == "files"}
    return rel, lab, file


def _feat_map(schema) -> dict:
    """{feature_name: Feature} for a schema — resolved once, reused for every row."""
    if schema is None:
        raise ValueError(
            "No schema on the record type — set up its features before syncing."
        )
    return {f.name: f for f in schema.members}


def _existing_by_ref(rec_type) -> dict:
    """{notion_uuid: ln.Record} for a type — ONE query, reused for upsert + write."""
    return {
        r.reference: r for r in ln.Record.filter(type=rec_type, reference_type="notion")
    }


def _resolved_map(uuids) -> dict:
    """{notion_uuid: ln.Record} for the UUIDs that resolve — ONE query, not N."""
    uuids = list(uuids)
    if not uuids:
        return {}
    return {
        r.reference: r
        for r in ln.Record.filter(reference__in=uuids, reference_type="notion")
    }


def _append_unique(values: list[str], value: str) -> None:
    if value not in values:
        values.append(value)


def _relation_target_from_feature(
    feature: Any,
    *,
    record_type_name: str,
    property_name: str,
) -> tuple[str, Any | None, str]:
    from lamindb.models.feature import parse_dtype

    dtype_str = getattr(feature, "_dtype_str", None)
    if not isinstance(dtype_str, str) or not dtype_str:
        raise ValueError(
            "Cannot sync relation values: feature "
            f"{record_type_name!r}.{property_name!r} has no dtype."
        )
    parsed = parse_dtype(dtype_str, check_exists=False)
    if len(parsed) != 1 or not parsed[0].get("list", False):
        raise ValueError(
            "Cannot sync relation values: feature "
            f"{record_type_name!r}.{property_name!r} must be typed as list[Record[type]] "
            f"or list[User]. "
            f"Got dtype={dtype_str!r}."
        )
    parsed_dtype = parsed[0]
    registry_str = parsed_dtype.get("registry_str")
    if registry_str == "User":
        return "user", None, "User"
    if registry_str != "Record":
        registry = parsed_dtype.get("registry")
        if registry is None:
            raise ValueError(
                "Cannot sync relation values: feature "
                f"{record_type_name!r}.{property_name!r} must target Record, User, or a Lamin registry type. "
                f"Got dtype={dtype_str!r}."
            )
        registry_name = getattr(registry, "__name__", str(registry))
        return "registry", registry, registry_name
    type_uid = parsed_dtype.get("type_uid")
    if not isinstance(type_uid, str) or not type_uid:
        raise ValueError(
            "Cannot sync relation values: feature "
            f"{record_type_name!r}.{property_name!r} must be a typed Record relation "
            "(e.g. list[Record[<type_uid>]]). Sync the source database metadata first."
        )
    qs = ln.Record.filter(uid=type_uid, is_type=True)
    count = qs.count()
    if count != 1:
        raise ValueError(
            "Cannot sync relation values: feature "
            f"{record_type_name!r}.{property_name!r} references type uid={type_uid!r} "
            f"but resolved {count} matching Record types."
        )
    target_type = qs.one()
    return "record", target_type, getattr(target_type, "name", "Record")


def _matches_target_type(record: Any, target_type: Any) -> bool:
    target_id = getattr(target_type, "id", None)
    if target_id is None:
        return True
    return getattr(record, "type_id", None) == target_id


def _relation_stub_name(reader: _NotionReader, notion_id: str) -> str:
    fallback = _compact_uuid(notion_id)
    try:
        payload = reader._call("GET", f"/pages/{_notion_api_id(notion_id)}")
    except Exception:
        return fallback
    name = _page_title(payload).strip()
    return name or fallback


def _relation_stub_detail(
    *,
    record_type_name: str,
    feature_name: str,
    target_type_name: str,
    stub_name: str,
    notion_id: str,
) -> str:
    return (
        f"{record_type_name} / {feature_name} -> {target_type_name}: "
        f"{stub_name} <- {_compact_uuid(notion_id)}"
    )


def _notion_user_name(reader: _NotionReader, notion_user_id: str) -> str | None:
    try:
        payload = reader._call("GET", f"/users/{_notion_api_id(notion_user_id)}")
    except Exception:
        return None
    name = payload.get("name")
    if isinstance(name, str) and name.strip():
        return name.strip()
    return None


def _notion_page_title(reader: _NotionReader, notion_page_id: str) -> str | None:
    try:
        payload = reader._call("GET", f"/pages/{_notion_api_id(notion_page_id)}")
    except Exception:
        return None
    title = _page_title(payload).strip()
    return title or None


def _notion_user_or_page_lookup(
    reader: _NotionReader, notion_id: str
) -> tuple[str | None, str]:
    cached_people_names = getattr(reader, "_notion_people_names_by_id", None)
    if isinstance(cached_people_names, dict):
        cached_name = cached_people_names.get(notion_id)
        if isinstance(cached_name, str) and cached_name.strip():
            return cached_name.strip(), "people_property_cache"

    notion_api_id = _notion_api_id(notion_id)

    user_reason = "no_user_name"
    try:
        user_payload = reader._call("GET", f"/users/{notion_api_id}")
        user_name = user_payload.get("name")
        if isinstance(user_name, str) and user_name.strip():
            return user_name.strip(), "user"
    except Exception as error:  # noqa: BLE001
        user_reason = type(error).__name__

    page_reason = "no_page_title"
    try:
        page_payload = reader._call("GET", f"/pages/{notion_api_id}")
        page_title = _page_title(page_payload).strip()
        if page_title:
            return page_title, "page"
    except Exception as error:  # noqa: BLE001
        page_reason = type(error).__name__

    return None, f"user:{user_reason};page:{page_reason}"


def _notion_user_or_page_name(reader: _NotionReader, notion_id: str) -> str | None:
    name, _ = _notion_user_or_page_lookup(reader, notion_id)
    return name


def _resolved_users_by_notion_id(
    reader: _NotionReader, notion_user_ids: list[str]
) -> dict[str, Any]:
    names_by_id: dict[str, str] = {}
    cached_people_names = getattr(reader, "_notion_people_names_by_id", None)
    for notion_user_id in notion_user_ids:
        if isinstance(cached_people_names, dict):
            cached_name = cached_people_names.get(notion_user_id)
            if isinstance(cached_name, str) and cached_name.strip():
                names_by_id[notion_user_id] = cached_name.strip()
                continue
        name = _notion_user_or_page_name(reader, notion_user_id)
        if name is not None:
            names_by_id[notion_user_id] = name
    if not names_by_id:
        return {}

    users_by_lower_name: dict[str, Any] = {}
    for name in sorted(set(names_by_id.values())):
        normalized_name = name.strip()
        normalized_key = normalized_name.lower()
        if not normalized_name:
            continue

        qs = ln.User.filter(name__iexact=normalized_name)
        name_count = qs.count()
        if name_count > 1:
            raise ValueError(
                f"Cannot sync relation values: ambiguous Lamin users for name {name!r}."
            )
        if name_count == 1:
            users_by_lower_name[normalized_key] = qs.one()

    return {
        notion_user_id: users_by_lower_name[name.lower()]
        for notion_user_id, name in names_by_id.items()
        if name.lower() in users_by_lower_name
    }


def _resolved_registry_records_by_notion_id(
    reader: _NotionReader, notion_page_ids: list[str], registry: Any
) -> dict[str, Any]:
    names_by_id: dict[str, str] = {}
    for notion_page_id in notion_page_ids:
        title = _notion_page_title(reader, notion_page_id)
        if title is not None:
            names_by_id[notion_page_id] = title
    if not names_by_id:
        return {}

    field_name = getattr(registry, "_name_field", "name")
    records_by_lower_name: dict[str, Any] = {}
    for name in sorted(set(names_by_id.values())):
        qs = registry.filter(**{f"{field_name}__iexact": name})
        count = qs.count()
        if count > 1:
            raise ValueError(
                "Cannot sync relation values: ambiguous Lamin records for "
                f"{registry.__name__}.{field_name}={name!r}."
            )
        if count == 1:
            records_by_lower_name[name.lower()] = qs.one()

    return {
        notion_page_id: records_by_lower_name[name.lower()]
        for notion_page_id, name in names_by_id.items()
        if name.lower() in records_by_lower_name
    }


def _registry_supports_relation_stubs(registry: Any) -> bool:
    registry_name = getattr(registry, "__name__", "")
    return registry_name in {"Project", "Reference"}


def _registry_stub_kwargs(
    registry: Any, field_name: str, stub_name: str, notion_id: str
) -> dict[str, str]:
    kwargs: dict[str, str] = {field_name: stub_name}
    if getattr(registry, "__name__", "") == "Project":
        compact_notion_id = _normalize_notion_id(notion_id) or notion_id
        kwargs["url"] = f"https://notion.so/laminlabs/{compact_notion_id}"
    return kwargs


def _resolve_relation_records_for_rows(
    reader: _NotionReader,
    rows: list[dict[str, Any]],
    rel: set[str],
    feat: dict[str, Any],
    prop_map: dict[str, str] | None,
    *,
    rec_type: Any,
    apply: bool,
    report: SyncReport | None,
) -> tuple[dict[str, Any], int]:
    prop_map = prop_map or {}
    relation_ids_by_prop: dict[str, set[str]] = {}
    for row in rows:
        for prop in rel:
            value = row.get(prop)
            if isinstance(value, list):
                ids = [
                    _normalize_notion_id(item) or item
                    for item in value
                    if isinstance(item, str) and item
                ]
            elif isinstance(value, str) and value:
                ids = [_normalize_notion_id(value) or value]
            else:
                ids = []
            if not ids:
                continue
            relation_ids_by_prop.setdefault(prop, set()).update(ids)

    resolved: dict[str, Any] = {}

    pending = 0
    record_type_name = getattr(rec_type, "name", str(rec_type))
    for prop in sorted(relation_ids_by_prop):
        feature_name = prop_map.get(prop, prop)
        feature = feat.get(feature_name)
        if feature is None:
            raise ValueError(
                f"Cannot sync relation values: missing feature mapping for property {prop!r}."
            )
        relation_target_kind, target_type, target_type_name = (
            _relation_target_from_feature(
                feature,
                record_type_name=record_type_name,
                property_name=feature_name,
            )
        )
        notion_ids = sorted(relation_ids_by_prop[prop])
        if relation_target_kind == "user":
            resolved_now = _resolved_users_by_notion_id(reader, notion_ids)
            for notion_id, user in resolved_now.items():
                resolved[notion_id] = user
            missing = [
                notion_id for notion_id in notion_ids if notion_id not in resolved_now
            ]
            existing_count = len(resolved_now)
            stub_create = 0
            unresolved_after = missing
        elif relation_target_kind == "registry":
            resolved_now = _resolved_registry_records_by_notion_id(
                reader, notion_ids, target_type
            )
            for notion_id, registry_record in resolved_now.items():
                resolved[notion_id] = registry_record
            missing = [
                notion_id for notion_id in notion_ids if notion_id not in resolved_now
            ]
            existing_count = len(resolved_now)
            stub_create = 0
            if target_type is not None and _registry_supports_relation_stubs(
                target_type
            ):
                field_name = getattr(target_type, "_name_field", "name")
                if apply:
                    for notion_id in missing:
                        stub_name = _relation_stub_name(reader, notion_id)
                        stub = target_type(
                            **_registry_stub_kwargs(
                                target_type,
                                field_name,
                                stub_name,
                                notion_id,
                            )
                        ).save()
                        resolved[notion_id] = stub
                        stub_create += 1
                        if report is not None:
                            _append_unique(
                                report.created_relation_stubs,
                                _relation_stub_detail(
                                    record_type_name=record_type_name,
                                    feature_name=feature_name,
                                    target_type_name=target_type_name,
                                    stub_name=stub_name,
                                    notion_id=notion_id,
                                ),
                            )
                else:
                    stub_create = len(missing)
                    if report is not None:
                        for notion_id in missing:
                            stub_name = _relation_stub_name(reader, notion_id)
                            _append_unique(
                                report.create_relation_stubs,
                                _relation_stub_detail(
                                    record_type_name=record_type_name,
                                    feature_name=feature_name,
                                    target_type_name=target_type_name,
                                    stub_name=stub_name,
                                    notion_id=notion_id,
                                ),
                            )
                unresolved_after = (
                    [notion_id for notion_id in notion_ids if notion_id not in resolved]
                    if apply
                    else []
                )
            else:
                unresolved_after = missing
        else:
            resolved_now = _resolved_map(notion_ids)
            for notion_id, record in resolved_now.items():
                resolved[notion_id] = record
            existing_count = len(resolved_now)
            missing = [
                notion_id for notion_id in notion_ids if notion_id not in resolved_now
            ]

            for notion_id, existing_record in resolved_now.items():
                if target_type is None:
                    continue
                if not _matches_target_type(existing_record, target_type):
                    raise ValueError(
                        "Cannot sync relation values: relation id "
                        f"{_compact_uuid(notion_id)!r} for {record_type_name!r}.{prop!r} "
                        f"resolves to record type {getattr(getattr(existing_record, 'type', None), 'name', None)!r} "
                        f"but expected {target_type_name!r}."
                    )

            stub_create = 0
            if missing:
                if apply:
                    for notion_id in missing:
                        stub_name = _relation_stub_name(reader, notion_id)
                        stub = ln.Record(
                            name=stub_name,
                            type=target_type,
                            reference=_normalize_notion_id(notion_id) or notion_id,
                            reference_type="notion",
                        ).save()
                        resolved[notion_id] = stub
                        stub_create += 1
                        if report is not None:
                            _append_unique(
                                report.created_relation_stubs,
                                _relation_stub_detail(
                                    record_type_name=record_type_name,
                                    feature_name=feature_name,
                                    target_type_name=target_type_name,
                                    stub_name=stub_name,
                                    notion_id=notion_id,
                                ),
                            )
                else:
                    stub_create = len(missing)
                    if report is not None:
                        for notion_id in missing:
                            stub_name = _relation_stub_name(reader, notion_id)
                            _append_unique(
                                report.create_relation_stubs,
                                _relation_stub_detail(
                                    record_type_name=record_type_name,
                                    feature_name=feature_name,
                                    target_type_name=target_type_name,
                                    stub_name=stub_name,
                                    notion_id=notion_id,
                                ),
                            )
            unresolved_after = (
                [notion_id for notion_id in notion_ids if notion_id not in resolved]
                if apply
                else []
            )
        pending += len(unresolved_after)
        if report is not None:
            _append_unique(
                report.relation_value_links,
                f"{record_type_name} / {prop}: "
                f"resolved_existing={existing_count}, "
                f"stub_create={stub_create}, "
                f"pending_unresolved={len(unresolved_after)}",
            )

    return resolved, pending


def _ensure_labels(names) -> None:
    """Create only the ULabels that don't already exist (exact-name match).

    Deterministic dedup by name: query the names already present in one shot,
    create just the missing ones. Never relies on the fuzzy similar-name search,
    so it can't blind-create a second 'intern' the way from_values(create=True)
    does when that search is disabled. (Pre-existing duplicates still need a
    one-time manual cleanup — this only stops NEW ones being made.)
    """
    names = {n for n in names if n}
    if not names:
        return
    existing = {u.name for u in ln.ULabel.filter(name__in=list(names))}
    for name in sorted(names - existing):
        ln.ULabel(name=name).save()


def _batch_labels(rows: list[dict], lab: set) -> None:
    """Ensure every distinct ULabel across all rows exists — one query, missing-only creates."""
    names = set()
    for row in rows:
        for p in lab:
            val = row.get(p)
            for n in val if isinstance(val, list) else [val]:
                if n:
                    names.add(n)
    _ensure_labels(names)


def _iter_file_urls(rows: list[dict], file_props: set[str]):
    for row in rows:
        notion_id = row.get("notion_id")
        for prop in file_props:
            value = row.get(prop)
            if isinstance(value, list):
                urls = value
            else:
                urls = [value]
            for url in urls:
                if isinstance(url, str) and url:
                    yield str(notion_id), prop, url


def _short_file_source(url: str, max_name_len: int = 48) -> str:
    parsed = urlparse(url)
    filename = Path(parsed.path).name or "file"
    if len(filename) > max_name_len:
        filename = filename[: max_name_len - 3] + "..."
    return filename


def _artifact_key_from_url(url: str) -> str:
    parsed = urlparse(url)
    filename = Path(parsed.path).name or "file"
    return f"notion_sync/{filename}"


def _iter_embedded_file_urls(content: str) -> set[str]:
    if not content:
        return set()
    urls: set[str] = set()
    urls.update(
        match.group(2) for match in MARKDOWN_IMAGE_LINK_PATTERN.finditer(content)
    )
    urls.update(match.group(2) for match in HTML_IMAGE_SRC_PATTERN.finditer(content))
    return {url for url in urls if isinstance(url, str) and url}


def _planned_embedded_file_transfers(
    markdown_content: str, notion_id: str
) -> tuple[dict[str, str], list[str]]:
    by_url: dict[str, str] = {}
    details: list[str] = []
    for url in sorted(_iter_embedded_file_urls(markdown_content)):
        detail = (
            f"{_short_file_source(url)} <- {_compact_uuid(notion_id)}:notes "
            '(key=None, kind="__easset__")'
        )
        details.append(detail)
        by_url.setdefault(url, detail)
    return by_url, details


def _storage_src_for_artifact(artifact: Any) -> str | None:
    protocol = None
    root_without_scheme = None
    root = getattr(getattr(ln.settings, "storage", None), "root", None)
    if root is not None:
        protocol = getattr(root, "protocol", None)
        root_str = str(root).strip()
        if protocol and root_str.startswith(f"{protocol}://"):
            root_without_scheme = root_str[len(protocol) + 3 :].rstrip("/")
        elif root_str:
            root_without_scheme = root_str.rstrip("/")
    artifact_path = getattr(artifact, "path", None)
    if artifact_path is None:
        return None
    path_str = str(artifact_path).strip()
    if not path_str:
        return None
    parsed = urlparse(path_str)
    artifact_protocol = parsed.scheme or protocol
    if not artifact_protocol:
        return None
    artifact_without_scheme = parsed.netloc + parsed.path
    artifact_without_scheme = artifact_without_scheme.strip("/")
    if not artifact_without_scheme:
        return None
    if (
        artifact_protocol == "s3"
        and root_without_scheme
        and artifact_without_scheme.startswith(root_without_scheme)
    ):
        rel = artifact_without_scheme[len(root_without_scheme) :]
        return f"/storage/s3/{root_without_scheme}%2F{rel}"
    if root_without_scheme and artifact_without_scheme.startswith(root_without_scheme):
        rel = artifact_without_scheme[len(root_without_scheme) :]
        encoded_rel = quote(rel, safe="/")
        return f"/storage/{artifact_protocol}/{root_without_scheme}{encoded_rel}"
    return f"/storage/{artifact_protocol}/{artifact_without_scheme}"


def _rewrite_embedded_file_refs(
    markdown_content: str, artifacts_by_url: dict[str, Any]
) -> str:
    if not markdown_content:
        return markdown_content

    def resolve(url: str) -> str:
        artifact = artifacts_by_url.get(url)
        if artifact is None:
            return url
        storage_src = _storage_src_for_artifact(artifact)
        return storage_src or url

    def replace_image(match: re.Match[str]) -> str:
        src = resolve(match.group(2))
        return f'\n\n<img width="{NOTION_EMBEDDED_IMAGE_WIDTH_PX}" src="{src}" />\n\n'

    def replace_html_image(match: re.Match[str]) -> str:
        return f"{match.group(1)}{resolve(match.group(2))}{match.group(3)}"

    content = MARKDOWN_IMAGE_LINK_PATTERN.sub(replace_image, markdown_content)
    content = HTML_IMAGE_SRC_PATTERN.sub(replace_html_image, content)
    return _normalize_markdown_block_spacing(content)


def _normalize_markdown_block_spacing(markdown_content: str) -> str:
    if not markdown_content:
        return markdown_content

    def classify(line: str, *, in_code_block: bool) -> str:
        stripped = line.strip()
        if stripped == "":
            return "blank"
        if re.match(r"^\s*```", line):
            return "code_fence"
        if in_code_block:
            return "code"
        if re.match(r"^\s{0,3}#{1,6}\s+", line):
            return "heading"
        if re.match(r"^\s*(?:[-+*]\s+|\d+\.\s+)", line):
            return "list"
        if re.match(r"^\s*(?:---|\*\*\*|___)\s*$", line):
            return "divider"
        if re.match(r"^\s*<img\b[^>]*>\s*$", line):
            return "html_img"
        return "paragraph"

    out: list[str] = []
    in_code_block = False
    prev_nonblank_kind: str | None = None

    for line in markdown_content.splitlines():
        kind = classify(line, in_code_block=in_code_block)

        needs_separator = False
        if kind != "blank":
            if prev_nonblank_kind in {"list"} and kind in {
                "heading",
                "paragraph",
                "html_img",
                "divider",
            }:
                needs_separator = True
            if kind in {"heading", "html_img", "divider"} and prev_nonblank_kind in {
                "paragraph",
                "list",
                "heading",
                "html_img",
                "divider",
            }:
                needs_separator = True
            if prev_nonblank_kind in {"heading", "html_img", "divider"} and kind in {
                "paragraph",
                "list",
            }:
                needs_separator = True

        if needs_separator and out and out[-1].strip() != "":
            out.append("")

        if kind == "blank":
            if out and out[-1].strip() != "":
                out.append("")
            continue

        out.append(line.rstrip())
        if kind == "code_fence":
            in_code_block = not in_code_block
        prev_nonblank_kind = kind if kind != "code_fence" else "code"

    return "\n".join(out).strip()


def _download_file_to_temp_path(url: str) -> str:
    parsed = urlparse(url)
    suffix = Path(parsed.path).suffix
    with tempfile.NamedTemporaryFile(
        delete=False,
        prefix="notion-sync-",
        suffix=suffix,
    ) as tmp:
        tmp_path = tmp.name
    with httpx.stream("GET", url, timeout=60, follow_redirects=True) as response:
        response.raise_for_status()
        with open(tmp_path, "wb") as handle:
            for chunk in response.iter_bytes():
                if chunk:
                    handle.write(chunk)
    return tmp_path


def _ensure_artifacts(
    file_urls: set[str],
    *,
    transfer_details_by_url: dict[str, str] | None = None,
    report: SyncReport | None = None,
    with_key: bool = True,
    kind: str | None = None,
    description: str | None = None,
) -> dict[str, Any]:
    """Create artifacts for Notion file URLs by downloading first."""
    urls = sorted({url for url in file_urls if isinstance(url, str) and url})
    if not urls:
        return {}
    artifacts: dict[str, Any] = {}
    for url in urls:
        tmp_path = None
        try:
            tmp_path = _download_file_to_temp_path(url)
            artifact_kwargs: dict[str, Any] = {}
            if with_key:
                artifact_kwargs["key"] = _artifact_key_from_url(url)
            if kind is not None:
                artifact_kwargs["kind"] = kind
            if description is not None:
                artifact_kwargs["description"] = description
            artifacts[url] = ln.Artifact(tmp_path, **artifact_kwargs).save()
            if report is not None and transfer_details_by_url is not None:
                detail = transfer_details_by_url.get(url, _short_file_source(url))
                if detail not in report.created_artifacts:
                    report.created_artifacts.append(detail)
        except (
            Exception
        ) as error:  # pragma: no cover - defensive against network/filesystem issues
            logger.warning(
                f"Could not register Notion file URL as Artifact: {url} ({error})"
            )
        finally:
            if tmp_path is not None:
                try:
                    Path(tmp_path).unlink()
                except OSError:
                    pass
    return artifacts


def _batch_artifacts(
    rows: list[dict],
    file_props: set[str],
    *,
    transfer_details_by_url: dict[str, str] | None = None,
    report: SyncReport | None = None,
    with_key: bool = True,
    kind: str | None = None,
    description: str | None = None,
) -> dict[str, Any]:
    urls = {url for _, _, url in _iter_file_urls(rows, file_props)}
    return _ensure_artifacts(
        urls,
        transfer_details_by_url=transfer_details_by_url,
        report=report,
        with_key=with_key,
        kind=kind,
        description=description,
    )


def _planned_missing_file_transfers(
    rows: list[dict], file_props: set[str]
) -> tuple[dict[str, str], list[str]]:
    by_url: dict[str, str] = {}
    details: list[str] = []
    for notion_id, prop, url in _iter_file_urls(rows, file_props):
        detail = f"{_short_file_source(url)} <- {_compact_uuid(notion_id)}:{prop}"
        details.append(detail)
        by_url.setdefault(url, detail)
    return by_url, details


def _planned_missing_embedded_transfers(
    reader: _NotionReader, rows: list[dict]
) -> tuple[dict[str, str], list[str]]:
    by_url: dict[str, str] = {}
    details: list[str] = []
    for row in rows:
        notion_id = row.get("notion_id")
        if not isinstance(notion_id, str) or not notion_id:
            continue
        try:
            markdown_content = reader.page_markdown(notion_id)
        except (
            Exception
        ) as error:  # pragma: no cover - defensive for API/network issues
            logger.warning(
                f"Could not inspect Notion page notes for embedded files: {notion_id} ({error})"
            )
            continue
        transfer_map, transfer_details = _planned_embedded_file_transfers(
            markdown_content, notion_id
        )
        for url, detail in transfer_map.items():
            by_url.setdefault(url, detail)
        details.extend(transfer_details)
    return by_url, details


def _row_values(
    row,
    rel,
    lab,
    file_props,
    feat,
    resolved,
    artifacts_by_url,
    prop_map,
    create_labels,
    skip_props,
):
    """Build the full {Feature: value} dict for one row. Returns (values, pending)."""
    prop_map = prop_map or {}
    values: dict[Any, Any] = {}
    pending = 0
    for prop, val in row.items():
        if (
            prop
            in (
                "notion_id",
                "created_time",
                "last_edited_time",
                "__notion_emoji__",
            )
            or prop in skip_props
            or val
            in (
                None,
                [],
                "",
            )
        ):
            continue
        f = feat.get(prop_map.get(prop, prop))
        if f is None:
            continue
        if prop in rel:
            raw_uuids = val if isinstance(val, list) else [val]
            uuids = [
                _normalize_notion_id(uuid_value) or uuid_value
                for uuid_value in raw_uuids
                if isinstance(uuid_value, str) and uuid_value
            ]
            hits = [resolved[u] for u in uuids if u in resolved]
            pending += sum(u not in resolved for u in uuids)
            if hits:
                values[f] = hits
        elif prop in lab:
            names = [n for n in (val if isinstance(val, list) else [val]) if n]
            if names:
                if (
                    create_labels
                ):  # standalone path; bulk path pre-creates via _batch_labels
                    _ensure_labels(names)
                values[f] = names if isinstance(val, list) else names[0]
        elif prop in file_props:
            urls = [u for u in (val if isinstance(val, list) else [val]) if u]
            hits = [artifacts_by_url[url] for url in urls if url in artifacts_by_url]
            if hits:
                values[f] = hits if isinstance(val, list) else hits[0]
        else:
            values[f] = _normalize_feature_value(f, val)
    return values, pending


def _normalize_feature_value(feature: Any, value: Any) -> Any:
    dtype = getattr(feature, "_dtype_str", None)
    if not isinstance(dtype, str):
        return value
    if dtype == "date" and isinstance(value, str):
        parsed_ts = _parse_notion_timestamp(value)
        if parsed_ts is not None:
            return parsed_ts.date()
        try:
            return date.fromisoformat(value)
        except ValueError:
            return value
    if dtype in {"datetime", "datetime64[ns, UTC]"} and isinstance(value, str):
        parsed_ts = _parse_notion_timestamp(value)
        if parsed_ts is not None:
            return parsed_ts
    return value


def _write(
    reader,
    rows,
    rec_type,
    spec,
    prop_map=None,
    by_id=None,
    *,
    transfer_details_by_url: dict[str, str] | None = None,
    report: SyncReport | None = None,
) -> dict:
    """Materialize every row of one database. Schema, kinds and labels resolved once."""
    rel, lab, file_props = _kinds(spec)
    _ensure_feature_itype_on_record_schema(rec_type)
    feat = _feat_map(rec_type.schema)
    internal_property_types = {
        "created_time",
        "last_edited_time",
        "created_by",
        "last_edited_by",
    }
    skip_props = {
        prop
        for prop, metadata in spec.items()
        if metadata["type"] in internal_property_types
    }

    _batch_labels(rows, lab)  # every ULabel created in one call
    artifacts_by_url = _batch_artifacts(
        rows,
        file_props,
        transfer_details_by_url=transfer_details_by_url,
        report=report,
    )
    markdown_by_notion_id: dict[str, str] = {}
    embedded_transfer_details_by_url: dict[str, str] = {}
    embedded_file_urls: set[str] = set()
    for row in rows:
        notion_id = row.get("notion_id")
        if not isinstance(notion_id, str) or not notion_id:
            continue
        markdown_content = reader.page_markdown(notion_id)
        markdown_by_notion_id[notion_id] = markdown_content
        per_page_transfer_map, _ = _planned_embedded_file_transfers(
            markdown_content, notion_id
        )
        for url, detail in per_page_transfer_map.items():
            embedded_file_urls.add(url)
            embedded_transfer_details_by_url.setdefault(url, detail)
    missing_embedded_urls = embedded_file_urls - set(artifacts_by_url)
    if missing_embedded_urls:
        artifacts_by_url.update(
            _ensure_artifacts(
                missing_embedded_urls,
                transfer_details_by_url=embedded_transfer_details_by_url,
                report=report,
                with_key=False,
                kind="__easset__",
                description="imported from Notion",
            )
        )

    relation_pending = 0
    if rel:
        resolved, relation_pending = _resolve_relation_records_for_rows(
            reader,
            rows,
            rel,
            feat,
            prop_map,
            rec_type=rec_type,
            apply=True,
            report=report,
        )
    else:
        resolved = {}

    if by_id is None:
        by_id = _existing_by_ref(rec_type)

    records = 0
    for row in rows:
        rec = by_id.get(row["notion_id"])
        if rec is None:  # not imported yet — nothing to write
            continue
        values, _ = _row_values(
            row,
            rel,
            lab,
            file_props,
            feat,
            resolved,
            artifacts_by_url,
            prop_map,
            create_labels=False,
            skip_props=skip_props,
        )
        rec.features.set_values(values)
        notion_id = row.get("notion_id")
        if isinstance(notion_id, str) and notion_id:
            markdown_content = markdown_by_notion_id.get(notion_id)
            if markdown_content is None:
                markdown_content = reader.page_markdown(notion_id)
            _attach_page_markdown(
                rec,
                _rewrite_embedded_file_refs(markdown_content, artifacts_by_url),
            )
        records += 1
    return {"records": records, "pending": relation_pending}


def _ensure_feature_itype_on_record_schema(rec_type: Any) -> None:
    schema = getattr(rec_type, "schema", None)
    if schema is None:
        return
    if getattr(schema, "itype", None) is not None:
        return
    members = getattr(schema, "members", None)
    if members is None:
        return
    first_member = None
    if hasattr(members, "all"):
        first_member = members.all().first()
    elif isinstance(members, list) and members:
        first_member = members[0]
    if first_member is None:
        return
    if first_member.__class__.__name__ != "Feature":
        return
    schema.itype = "Feature"
    save_fn = getattr(schema, "save", None)
    if callable(save_fn):
        save_fn(update_fields=["itype"])
        logger.warning(
            "notion sync metadata: repaired schema itype to 'Feature' "
            f"for record type {getattr(rec_type, 'name', rec_type)!r}"
        )


def _attach_page_markdown(record: Any, markdown_content: str) -> None:
    content = markdown_content.strip()
    if not content:
        return
    if getattr(record, "notes", None) == content:
        return
    recordblock = ln.models.RecordBlock(
        record=record,
        content=content,
        kind="readme",
    ).save()
    record.ablocks.add(recordblock, bulk=False)


def _attach_project_markdown(project: Any, markdown_content: str) -> None:
    content = markdown_content.strip()
    if not content:
        return
    if getattr(project, "notes", None) == content:
        return
    projectblock = ln.models.ProjectBlock(
        project=project,
        content=content,
        kind="readme",
    ).save()
    project.ablocks.add(projectblock, bulk=False)


def _upsert_all(rec_type, rows) -> dict:
    """Upsert every row of one database against a single existing-record map.

    Returns the {notion_uuid: record} map (existing + newly created), ready to
    hand to :func:`_write` so it never re-queries.
    """
    by_id = _existing_by_ref(rec_type)  # ONE query, not one per row
    for row in rows:
        raw_nid, name = row["notion_id"], row.get("name")
        if not isinstance(raw_nid, str) or not raw_nid:
            continue
        nid = _normalize_notion_id(raw_nid) or raw_nid
        row_emoji = row.get("__notion_emoji__")
        created_at = _parse_notion_timestamp(row.get("created_time"))
        updated_at = _parse_notion_timestamp(row.get("last_edited_time"))
        if created_at is None:
            created_at = updated_at
        if updated_at is None:
            updated_at = created_at
        rec = by_id.get(nid)
        if rec is None:
            rec = ln.Record(
                name=name or None,
                type=rec_type,
                reference=nid,
                reference_type="notion",
            )
            rec._aux = _NotionSyncer._merge_aux_with_emoji(None, row_emoji)
            if created_at is not None:
                rec.created_at = created_at
            if updated_at is not None:
                rec.updated_at = updated_at
            by_id[nid] = rec.save()
            continue

        changed_fields: list[str] = []
        if getattr(rec, "reference", None) != nid:
            rec.reference = nid
            changed_fields.append("reference")
        if name and rec.name != name:
            rec.name = name
            changed_fields.append("name")
        if created_at is not None and _normalized_timestamp(
            rec.created_at
        ) != _normalized_timestamp(created_at):
            rec.created_at = created_at
            changed_fields.append("created_at")
        if updated_at is not None and _normalized_timestamp(
            rec.updated_at
        ) != _normalized_timestamp(updated_at):
            rec.updated_at = updated_at
            changed_fields.append("updated_at")
        merged_aux = _NotionSyncer._merge_aux_with_emoji(
            getattr(rec, "_aux", None), row_emoji
        )
        if getattr(rec, "_aux", None) != merged_aux:
            rec._aux = merged_aux
            changed_fields.append("_aux")
        if changed_fields:
            rec.save(update_fields=changed_fields)
    return by_id


class _NotionSyncer:
    """Sync Notion page trees to typed LaminDB records."""

    def __init__(self, token: str | None = None) -> None:
        token = token or os.getenv("NOTION_TOKEN")
        if not token:
            raise ValueError("Pass token=... or set NOTION_TOKEN.")
        self.reader = _NotionReader(token=token)
        self._parent_page_emojis: dict[str, str | None] = {}
        self._database_parent_pages: dict[str, str] = {}
        self._seed_page_ids_by_database: dict[str, set[str]] = {}
        self._parent_page_parents: dict[str, str] = {}
        self._formula_dtype_cache: dict[tuple[str, str], tuple[str, Any]] = {}

    @staticmethod
    def _parse_formula_dtype_choice(choice: str) -> tuple[str, Any] | None:
        normalized = choice.strip().lower()
        mapping: dict[str, tuple[str, Any]] = {
            "bool": ("bool", bool),
            "boolean": ("bool", bool),
            "num": ("num", "num"),
            "number": ("num", "num"),
            "str": ("str", str),
            "string": ("str", str),
            "url": ("url", "url"),
            "datetime": ("datetime64[ns, UTC]", datetime),
            "datetime64[ns, utc]": ("datetime64[ns, UTC]", datetime),
            "date": ("datetime64[ns, UTC]", datetime),
        }
        return mapping.get(normalized)

    def _existing_formula_dtype_choice(
        self, db_name: str, property_name: str
    ) -> tuple[str, Any] | None:
        schema_qs = ln.Schema.filter(name=db_name)
        if schema_qs.count() != 1:
            return None
        schema = schema_qs.one_or_none()
        if schema is None:
            return None
        members = schema.members
        if hasattr(members, "filter"):
            candidates = list(members.filter(name__iexact=property_name))
        else:
            candidates = [
                feature
                for feature in members
                if isinstance(getattr(feature, "name", None), str)
                and feature.name.lower() == property_name.lower()
            ]
        existing_feature = self._pick_unique(candidates)
        if existing_feature is None:
            return None
        existing_dtype = getattr(existing_feature, "dtype_as_str", None) or getattr(
            existing_feature, "_dtype_str", None
        )
        if not isinstance(existing_dtype, str):
            return None
        return self._parse_formula_dtype_choice(existing_dtype)

    def _dtype_from_formula_property(
        self,
        db_name: str,
        property_name: str,
        property_spec: dict[str, Any],
    ) -> tuple[str, Any]:
        cache_key = (db_name, property_name)
        cached = self._formula_dtype_cache.get(cache_key)
        if cached is not None:
            return cached
        existing = self._existing_formula_dtype_choice(db_name, property_name)
        if existing is not None:
            logger.important(
                "formula dtype: reusing existing schema dtype "
                f"db={db_name!r}, property={property_name!r}, dtype={existing[0]!r}"
            )
            self._formula_dtype_cache[cache_key] = existing
            return existing
        expression = property_spec.get("formula_expression")
        if not isinstance(expression, str) or not expression:
            expression = "<formula expression unavailable from Notion API>"
        logger.important(
            f"notion sync formula dtype prompt: db={db_name!r}, property={property_name!r}"
        )
        RICH_CONSOLE.print(
            f"[bold yellow]notion formula[/] {db_name}.{property_name} = {expression}",
            markup=True,
            highlight=False,
        )
        if sys.stdin is None or not sys.stdin.isatty():
            logger.warning(
                f"notion sync formula dtype: non-interactive session, defaulting {db_name}.{property_name} to str"
            )
            resolved = ("str", str)
            self._formula_dtype_cache[cache_key] = resolved
            return resolved
        while True:
            selected = input(
                f"Choose Lamin dtype for {db_name}.{property_name} "
                "[bool/num/str/url/datetime] (default=str): "
            ).strip()
            if selected == "":
                resolved = ("str", str)
                self._formula_dtype_cache[cache_key] = resolved
                return resolved
            parsed = self._parse_formula_dtype_choice(selected)
            if parsed is not None:
                self._formula_dtype_cache[cache_key] = parsed
                return parsed
            RICH_CONSOLE.print(
                "[yellow]Invalid dtype. Use one of: bool, num, str, url, datetime.[/]",
                markup=True,
                highlight=False,
            )

    def _safe_call(self, path: str) -> dict | None:
        try:
            return self.reader._call("GET", path)
        except LookupError:
            return None
        except httpx.HTTPStatusError as error:
            # Parent probing can hit 400 when an ID is valid but not for the
            # probed endpoint (e.g. page ID on /databases/{id}); treat as miss.
            if error.response is not None and error.response.status_code == 400:
                return None
            raise

    def _iter_block_children(self, block_id: str) -> list[dict]:
        children: list[dict] = []
        cursor: str | None = None
        while True:
            params: dict[str, Any] = {"page_size": 100}
            if cursor:
                params["start_cursor"] = cursor
            payload = self.reader._call(
                "GET", f"/blocks/{_notion_api_id(block_id)}/children", params=params
            )
            children.extend(payload.get("results", []))
            if not payload.get("has_more"):
                return children
            cursor = payload.get("next_cursor")

    def _collect_databases_from_block(
        self,
        block_id: str,
        seen: set[str],
        out: set[str],
        *,
        parent_page_id: str | None = None,
        depth: int | None = None,
    ) -> None:
        """Collect child databases up to `depth` levels below this block.

        ``depth=None`` walks the whole subtree. ``depth<=0`` does not enter it.
        A child database on this block is one level; databases nested under it
        consume the remaining levels.
        """
        if depth is not None and depth <= 0:
            return
        next_depth = None if depth is None else depth - 1
        for block in self._iter_block_children(block_id):
            bid = block.get("id")
            if bid and bid in seen:
                continue
            if bid:
                seen.add(bid)
            if block.get("type") == "child_database" and bid:
                out.add(bid)
                normalized_db_id = _normalize_notion_id(bid)
                if parent_page_id is not None and normalized_db_id is not None:
                    self._database_parent_pages[normalized_db_id] = parent_page_id
            if block.get("has_children") and bid:
                self._collect_databases_from_block(
                    bid,
                    seen,
                    out,
                    parent_page_id=parent_page_id,
                    depth=next_depth,
                )

    def _collect_database_ids(
        self, parents: list[str], depth: int | None = None
    ) -> tuple[set[str], dict[str, str]]:
        if depth is not None and depth < 0:
            raise ValueError("depth must be >= 0 when provided.")
        database_ids: set[str] = set()
        parent_pages: dict[str, str] = {}
        parent_page_emojis: dict[str, str | None] = {}
        self._database_parent_pages = {}
        self._seed_page_ids_by_database = {}
        self._parent_page_parents = {}
        seen_blocks: set[str] = set()
        for parent in parents:
            db_payload = self._safe_call(f"/databases/{parent}")
            if db_payload is not None:
                database_ids.add(parent)
                normalized_db_id = _normalize_notion_id(parent) or parent
                parent_page_id = self._database_parent_page_id(db_payload)
                if parent_page_id is not None:
                    self._database_parent_pages[normalized_db_id] = parent_page_id
                    if parent_page_id not in parent_pages:
                        parent_page_payload = self._safe_call(
                            f"/pages/{parent_page_id}"
                        )
                        if parent_page_payload is not None:
                            parent_title = _page_title(
                                parent_page_payload
                            ).strip() or _compact_uuid(parent_page_id)
                            parent_pages[parent_page_id] = parent_title
                            parent_page_emojis[parent_page_id] = self._database_emoji(
                                parent_page_payload
                            )
                # recurse through rows as pages to discover nested child databases
                if depth == 0:
                    continue
                for row in self.reader.rows(parent):
                    notion_id = row.get("notion_id")
                    if notion_id:
                        self._collect_databases_from_block(
                            notion_id,
                            seen_blocks,
                            database_ids,
                            parent_page_id=None,
                            depth=depth,
                        )
                continue
            page_payload = self._safe_call(f"/pages/{parent}")
            if page_payload is None:
                raise LookupError(
                    f"Parent {_compact_uuid(parent)!r} is neither a readable database nor page."
                )
            parent_id = _normalize_notion_id(parent) or parent
            parent_database_id = self._page_parent_database_id(page_payload)
            if parent_database_id is not None:
                normalized_db_id = (
                    _normalize_notion_id(parent_database_id) or parent_database_id
                )
                database_ids.add(normalized_db_id)
                self._seed_page_ids_by_database.setdefault(normalized_db_id, set()).add(
                    parent_id
                )
                db_payload = self._safe_call(f"/databases/{parent_database_id}")
                if db_payload is not None:
                    db_parent_page_id = self._database_parent_page_id(db_payload)
                    if db_parent_page_id is not None:
                        self._database_parent_pages[normalized_db_id] = (
                            db_parent_page_id
                        )
                if depth == 0:
                    continue
                self._collect_databases_from_block(
                    parent,
                    seen_blocks,
                    database_ids,
                    parent_page_id=None,
                    depth=depth,
                )
                continue

            parent_title = _page_title(page_payload).strip() or _compact_uuid(parent)
            parent_pages[parent_id] = parent_title
            parent_page_emojis[parent_id] = self._database_emoji(page_payload)
            ancestor_page_id = self._database_parent_page_id(page_payload)
            if ancestor_page_id is not None:
                self._parent_page_parents[parent_id] = ancestor_page_id
            if ancestor_page_id is not None and ancestor_page_id not in parent_pages:
                ancestor_page_payload = self._safe_call(f"/pages/{ancestor_page_id}")
                if ancestor_page_payload is not None:
                    ancestor_title = _page_title(
                        ancestor_page_payload
                    ).strip() or _compact_uuid(ancestor_page_id)
                    parent_pages[ancestor_page_id] = ancestor_title
                    parent_page_emojis[ancestor_page_id] = self._database_emoji(
                        ancestor_page_payload
                    )
            if depth == 0:
                continue
            self._collect_databases_from_block(
                parent,
                seen_blocks,
                database_ids,
                parent_page_id=parent_id,
                depth=depth,
            )
        self._parent_page_emojis = parent_page_emojis
        return database_ids, parent_pages

    def _rows_for_seed_pages(
        self, database_id: str, page_ids: set[str], *, include_page_emoji: bool = False
    ) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        normalized_database_id = _normalize_notion_id(database_id) or database_id
        for page_id in sorted(page_ids):
            page_payload = self._safe_call(f"/pages/{_notion_api_id(page_id)}")
            if page_payload is None:
                continue
            page_database_id = self._page_parent_database_id(page_payload)
            if (page_database_id or "") != normalized_database_id:
                continue
            row: dict[str, Any] = {
                "notion_id": _normalize_notion_id(page_payload.get("id")),
                "created_time": page_payload.get("created_time"),
                "last_edited_time": page_payload.get("last_edited_time"),
            }
            properties = page_payload.get("properties", {})
            if isinstance(properties, dict):
                for name, prop in properties.items():
                    if isinstance(prop, dict):
                        row[name] = _flatten(
                            prop,
                            people_name_cache=self.reader._notion_people_names_by_id,
                        )
            if include_page_emoji:
                row["__notion_emoji__"] = _extract_emoji(page_payload)
            if row["notion_id"]:
                rows.append(row)
        return rows

    def _resolve_or_create_type_by_name(
        self,
        name: str,
        *,
        emoji: str | None = None,
        apply: bool,
        report: SyncReport,
    ):
        qs = ln.Record.filter(name=name, is_type=True)
        count = qs.count()
        if count == 0:
            if apply:
                record_kwargs: dict[str, Any] = {"name": name, "is_type": True}
                aux = self._merge_aux_with_emoji(None, emoji)
                if aux is not None:
                    record_kwargs["_aux"] = aux
                rec_type = ln.Record(**record_kwargs).save()
                if name not in report.created_record_types:
                    report.created_record_types.append(name)
                return rec_type
            elif name not in report.create_record_types:
                report.create_record_types.append(name)
            return None
        if count > 1:
            raise ValueError(
                f"Ambiguous Lamin record type name {name!r}: found {count} matches."
            )
        rec_type = qs.one()
        if apply and emoji is not None:
            merged_aux = self._merge_aux_with_emoji(
                getattr(rec_type, "_aux", None), emoji
            )
            if getattr(rec_type, "_aux", None) != merged_aux:
                rec_type._aux = merged_aux
                rec_type.save(update_fields=["_aux"])
        return rec_type

    @staticmethod
    def _database_parent_page_id(payload: dict) -> str | None:
        parent = payload.get("parent")
        if not isinstance(parent, dict):
            return None
        if parent.get("type") != "page_id":
            return None
        page_id = parent.get("page_id")
        return _normalize_notion_id(page_id)

    def _page_parent_database_id(self, payload: dict) -> str | None:
        parent = payload.get("parent")
        if not isinstance(parent, dict):
            return None
        parent_type = parent.get("type")
        if parent_type == "database_id":
            return _normalize_notion_id(parent.get("database_id"))
        if parent_type == "data_source_id":
            data_source_id = parent.get("data_source_id")
            if not isinstance(data_source_id, str) or not data_source_id.strip():
                return None
            data_source_payload = self._safe_call(f"/data_sources/{data_source_id}")
            if data_source_payload is None:
                return None
            ds_parent = data_source_payload.get("parent")
            if not isinstance(ds_parent, dict):
                return None
            if ds_parent.get("type") != "database_id":
                return None
            return _normalize_notion_id(ds_parent.get("database_id"))
        return None

    @staticmethod
    def _database_title(payload: dict, fallback: str) -> str:
        title = payload.get("title") or []
        text = "".join(part.get("plain_text", "") for part in title).strip()
        return text or fallback

    @staticmethod
    def _database_description(payload: dict) -> str | None:
        description = payload.get("description") or []
        if not isinstance(description, list):
            return None
        text = "".join(part.get("plain_text", "") for part in description).strip()
        return text or None

    @staticmethod
    def _database_emoji(payload: dict) -> str | None:
        return _extract_emoji(payload)

    @staticmethod
    def _merge_aux_with_emoji(aux: Any, emoji: str | None) -> dict[str, Any] | None:
        merged = dict(aux) if isinstance(aux, dict) else {}
        if emoji is None:
            merged.pop("ei", None)
        else:
            merged["ei"] = emoji
        return merged or None

    @staticmethod
    def _schema_feature_names(rec_type) -> set[str]:
        schema = rec_type.schema
        if schema is None:
            raise ValueError(
                f"Record type {rec_type.name!r} has no schema. Add a schema before syncing."
            )
        return {feature.name for feature in schema.members}

    @staticmethod
    def _record_type_schema_attach_detail(rec_type_name: str, schema_name: str) -> str:
        return f"{rec_type_name}: attach schema {schema_name}"

    @staticmethod
    def _feature_dtype_from_notion_type(notion_type: str):
        if notion_type == "number":
            return "num"
        if notion_type == "url":
            return "url"
        if notion_type == "checkbox":
            return bool
        if notion_type in {"created_time", "last_edited_time"}:
            return datetime
        if notion_type in {"created_by", "last_edited_by"}:
            return ln.User
        if notion_type == "people":
            return list[ln.User]
        if notion_type in {"multi_select", "relation"}:
            return list[str]
        if notion_type == "files":
            return list[ln.Artifact]
        return str

    @staticmethod
    def _feature_dtype_label_from_notion_type(notion_type: str) -> str:
        if notion_type == "number":
            return "num"
        if notion_type == "url":
            return "url"
        if notion_type == "checkbox":
            return "bool"
        if notion_type in {"created_time", "last_edited_time"}:
            return "datetime64[ns, UTC]"
        if notion_type in {"created_by", "last_edited_by"}:
            return "User"
        if notion_type == "people":
            return "list[User]"
        if notion_type in {"multi_select", "relation"}:
            return "list[str]"
        if notion_type == "files":
            return "list[Artifact]"
        return "str"

    @staticmethod
    def _name_candidates(name: str) -> list[str]:
        raw = (name or "").strip()
        if not raw:
            return []
        base = raw.replace("_", " ")
        variants = [raw, base, base.title()]
        if base.endswith("s"):
            singular = base[:-1].strip()
            if singular:
                variants.extend([singular, singular.title()])
        else:
            plural = f"{base}s"
            variants.extend([plural, plural.title()])
        # Preserve order while deduplicating.
        return list(dict.fromkeys(v for v in variants if v))

    @staticmethod
    def _strict_name_candidates(name: str) -> list[str]:
        raw = (name or "").strip()
        if not raw:
            return []
        base = raw.replace("_", " ")
        variants = [raw, base, base.title()]
        return list(dict.fromkeys(v for v in variants if v))

    @staticmethod
    def _pick_unique(records: list[Any]) -> Any | None:
        unique_by_id: dict[Any, Any] = {}
        for record in records:
            record_id = getattr(record, "id", None)
            key = record_id if record_id is not None else id(record)
            unique_by_id[key] = record
        if len(unique_by_id) == 1:
            return next(iter(unique_by_id.values()))
        return None

    def _resolve_record_type_by_name_candidates(self, names: list[str]) -> Any | None:
        matches: list[Any] = []
        for candidate in names:
            qs = ln.Record.filter(name__iexact=candidate, is_type=True)
            if qs.count() == 1:
                matches.append(qs.one())
        return self._pick_unique(matches)

    def _resolve_ulabel_type_by_name_candidates(
        self, names: list[str], *, parent_type: Any | None = None
    ) -> Any | None:
        matches: list[Any] = []
        for candidate in names:
            filters: dict[str, Any] = {"name__iexact": candidate, "is_type": True}
            if parent_type is None:
                filters["type__isnull"] = True
            else:
                filters["type"] = parent_type
            qs = ln.ULabel.filter(**filters)
            if qs.count() == 1:
                matches.append(qs.one())
        return self._pick_unique(matches)

    @staticmethod
    def _list_dtype_for(dynamic_type: Any) -> Any:
        # dynamic_type is resolved at runtime from DB records, so mypy cannot
        # validate it as a static type argument.
        return list[dynamic_type]  # type: ignore[valid-type]

    @staticmethod
    def _singularize(value: str) -> str:
        base = value.strip()
        if base.endswith("s") and len(base) > 1:
            return base[:-1]
        return base

    @staticmethod
    def _pluralize(value: str) -> str:
        base = value.strip()
        if not base:
            return base
        if base.endswith("y") and len(base) > 1 and base[-2].lower() not in "aeiou":
            return f"{base[:-1]}ies"
        if base.endswith("s"):
            return base
        return f"{base}s"

    def _label_type_name_candidates(
        self, db_name: str, property_name: str
    ) -> list[str]:
        db_candidates = self._name_candidates(db_name)
        prop_candidates = self._name_candidates(property_name)
        singular_db = self._singularize(db_name.replace("_", " "))
        singular_db_candidates = self._name_candidates(singular_db)
        combined: list[str] = []
        for db_candidate in [*db_candidates, *singular_db_candidates]:
            for prop_candidate in prop_candidates:
                combined.append(f"{db_candidate} {prop_candidate}".strip())
        return list(dict.fromkeys([*combined, *prop_candidates]))

    def _default_label_type_name(self, db_name: str, property_name: str) -> str:
        prop = property_name.replace(" ", "_").strip().lower()
        return self._pluralize(prop)

    def _resolve_or_plan_parent_ulabel_type(
        self,
        db_name: str,
        *,
        apply: bool,
        report: SyncReport | None,
    ) -> Any | None:
        parent = self._resolve_ulabel_type_by_name_candidates(
            [db_name], parent_type=None
        )
        if parent is not None:
            return parent
        if report is not None:
            key = "created_ulabel_types" if apply else "create_ulabel_types"
            self._append_unique(getattr(report, key), db_name)
        if apply:
            return ln.ULabel(name=db_name, is_type=True).save()
        return None

    def _resolve_or_plan_ulabel_type(
        self,
        db_name: str,
        property_name: str,
        *,
        apply: bool,
        report: SyncReport | None,
    ) -> tuple[Any | None, str, str]:
        parent_type = self._resolve_or_plan_parent_ulabel_type(
            db_name, apply=apply, report=report
        )
        default_name = self._default_label_type_name(db_name, property_name)
        candidates = [default_name, *self._name_candidates(property_name)]
        candidates = list(dict.fromkeys(c for c in candidates if c))
        label_type = self._resolve_ulabel_type_by_name_candidates(
            candidates, parent_type=parent_type
        )
        label_type_path = f"{db_name} / {default_name}"
        if label_type is not None:
            return label_type, label_type.name, label_type_path
        if report is not None:
            key = "created_ulabel_types" if apply else "create_ulabel_types"
            self._append_unique(getattr(report, key), label_type_path)
        if apply:
            created = ln.ULabel(
                name=default_name, is_type=True, type=parent_type
            ).save()
            return created, created.name, label_type_path
        return None, default_name, label_type_path

    def _plan_or_create_ulabels(
        self,
        label_type: Any | None,
        label_type_path: str,
        choices: list[str] | None,
        *,
        apply: bool,
        report: SyncReport | None,
    ) -> None:
        if report is None or not choices:
            return
        ordered_unique = sorted({value for value in choices if value})
        if not ordered_unique:
            return
        existing: set[str] = set()
        if label_type is not None:
            existing = set(
                ln.ULabel.filter(name__in=ordered_unique, type=label_type).values_list(
                    "name", flat=True
                )
            )
        missing = [value for value in ordered_unique if value not in existing]
        for value in missing:
            detail = f"{label_type_path} / {value}"
            key = "created_ulabels" if apply else "create_ulabels"
            self._append_unique(getattr(report, key), detail)
        if apply and label_type is not None:
            for value in missing:
                ln.ULabel(name=value, type=label_type).save()

    def _relation_target_name_candidates(self, target: str) -> list[str]:
        normalized_target = _normalize_notion_id(target) or target

        db_payload = self._safe_call(f"/databases/{target}")
        if db_payload is not None:
            db_name = self._database_title(db_payload, fallback=normalized_target)
            return self._strict_name_candidates(db_name)

        data_source_payload = self._safe_call(f"/data_sources/{target}")
        if data_source_payload is None:
            return []

        names: list[str] = []
        ds_name = data_source_payload.get("name")
        if isinstance(ds_name, str) and ds_name.strip():
            names.extend(self._strict_name_candidates(ds_name.strip()))

        parent = data_source_payload.get("parent")
        if isinstance(parent, dict):
            parent_db_id = parent.get("database_id")
            if isinstance(parent_db_id, str) and parent_db_id.strip():
                parent_db_payload = self._safe_call(f"/databases/{parent_db_id}")
                if parent_db_payload is not None:
                    db_name = self._database_title(
                        parent_db_payload, fallback=parent_db_id
                    )
                    names.extend(self._strict_name_candidates(db_name))

        return list(dict.fromkeys(names))

    @staticmethod
    def _preferred_record_type_name(
        names: list[str], *, prefer_plural: bool = True
    ) -> str | None:
        candidates = [
            name for name in dict.fromkeys(names) if isinstance(name, str) and name
        ]
        if not candidates:
            return None
        if not prefer_plural:
            titled = [name for name in candidates if name[:1].isupper()]
            if titled:
                return titled[0]
            return candidates[0]
        plural_titled = [
            name for name in candidates if name[:1].isupper() and name.endswith("s")
        ]
        if plural_titled:
            return plural_titled[0]
        titled = [name for name in candidates if name[:1].isupper()]
        if titled:
            return titled[0]
        plural = [name for name in candidates if name.endswith("s")]
        if plural:
            return plural[0]
        return candidates[0]

    def _resolve_or_plan_relation_record_type(
        self,
        names: list[str],
        *,
        apply: bool,
        report: SyncReport | None,
        prefer_plural: bool = True,
    ) -> tuple[Any | None, str | None]:
        deduped = list(dict.fromkeys(names))
        normalized = {
            name.strip().lower().replace("_", " ")
            for name in deduped
            if isinstance(name, str) and name.strip()
        }
        if normalized & {"project", "projects"}:
            return ln.Project, "Project"
        if normalized & {"reference", "references"}:
            return ln.Reference, "Reference"
        relation_type = self._resolve_record_type_by_name_candidates(deduped)
        if relation_type is not None:
            return relation_type, relation_type.name
        planned_name = self._preferred_record_type_name(
            deduped, prefer_plural=prefer_plural
        )
        if planned_name is None or report is None:
            return None, None
        created = self._resolve_or_create_type_by_name(
            planned_name, apply=apply, report=report
        )
        if created is not None:
            return created, created.name
        return None, planned_name

    def _dtype_from_notion_property(
        self,
        db_name: str,
        property_name: str,
        property_spec: dict[str, Any],
        *,
        apply: bool = False,
        report: SyncReport | None = None,
    ) -> tuple[str, Any]:
        notion_type = property_spec["type"]
        if notion_type == "formula":
            return self._dtype_from_formula_property(
                db_name, property_name, property_spec
            )
        if notion_type == "people":
            return "list[User]", list[ln.User]
        if notion_type in {"select", "status", "multi_select"}:
            label_type, label_type_name, label_type_path = (
                self._resolve_or_plan_ulabel_type(
                    db_name,
                    property_name,
                    apply=apply,
                    report=report,
                )
            )
            self._plan_or_create_ulabels(
                label_type,
                label_type_path,
                property_spec.get("choices"),
                apply=apply,
                report=report,
            )
            if notion_type == "multi_select":
                if label_type is not None:
                    return f"list[{label_type.name}]", self._list_dtype_for(label_type)
                return f"list[{label_type_name}]", list[ln.ULabel]
            if label_type is not None:
                return label_type.name, label_type
            return label_type_name, ln.ULabel
        if notion_type == "relation":
            target = property_spec.get("target")
            if isinstance(target, str) and target:
                target_names = self._relation_target_name_candidates(target)
                relation_type, planned_name = (
                    self._resolve_or_plan_relation_record_type(
                        target_names, apply=apply, report=report, prefer_plural=False
                    )
                )
                if relation_type is not None:
                    if relation_type in {ln.Project, ln.Reference}:
                        logger.important(
                            "notion sync metadata: "
                            f"{db_name}.{property_name} relation maps to "
                            f"LaminDB {relation_type.__name__} registry"
                        )
                    relation_type_name = getattr(relation_type, "name", None)
                    if not isinstance(relation_type_name, str):
                        relation_type_name = getattr(
                            relation_type, "__name__", str(relation_type)
                        )
                    return f"list[{relation_type_name}]", self._list_dtype_for(
                        relation_type
                    )
                if planned_name is not None:
                    return f"list[{planned_name}]", list[ln.Record]
                # Fallback to the local property name only (not dual synced-name)
                # to recover obvious mappings like presentations->Presentations
                # while avoiding cross-side self-type mis-inference.
                relation_type, planned_name = (
                    self._resolve_or_plan_relation_record_type(
                        self._name_candidates(property_name),
                        apply=apply,
                        report=report,
                    )
                )
                if relation_type is not None:
                    return f"list[{relation_type.name}]", self._list_dtype_for(
                        relation_type
                    )
                if planned_name is not None:
                    return f"list[{planned_name}]", list[ln.Record]
                return "list[str]", list[str]

            names = self._name_candidates(property_name)
            dual = property_spec.get("dual")
            if isinstance(dual, dict):
                synced_name = dual.get("synced_property_name")
                if isinstance(synced_name, str) and synced_name:
                    names.extend(self._name_candidates(synced_name))
            special_relation_type, _ = self._resolve_or_plan_relation_record_type(
                list(dict.fromkeys(names)),
                apply=apply,
                report=None,
            )
            if special_relation_type is not None and special_relation_type in {
                ln.Project,
                ln.Reference,
            }:
                logger.important(
                    "notion sync metadata: "
                    f"{db_name}.{property_name} relation maps to "
                    f"LaminDB {special_relation_type.__name__} registry"
                )
                return (
                    f"list[{special_relation_type.__name__}]",
                    self._list_dtype_for(special_relation_type),
                )
            relation_type = self._resolve_record_type_by_name_candidates(
                list(dict.fromkeys(names))
            )
            if relation_type is not None:
                relation_type_name = getattr(relation_type, "name", None)
                if not isinstance(relation_type_name, str):
                    relation_type_name = getattr(
                        relation_type, "__name__", str(relation_type)
                    )
                return f"list[{relation_type_name}]", self._list_dtype_for(
                    relation_type
                )
            return "list[str]", list[str]
        return (
            self._feature_dtype_label_from_notion_type(notion_type),
            self._feature_dtype_from_notion_type(notion_type),
        )

    @staticmethod
    def _record_field_mapping_for_notion_type(notion_type: str) -> str | None:
        if notion_type == "title":
            return "name"
        if notion_type == "created_time":
            return "created_at"
        if notion_type == "last_edited_time":
            return "updated_at"
        if notion_type == "created_by":
            return "created_by"
        return None

    @staticmethod
    def _record_field_mapping_for_property_name(property_name: str) -> str | None:
        normalized = property_name.strip().lower().replace(" ", "_")
        if normalized in {"summary", "description"}:
            return "description"
        return None

    def _record_field_mappings_from_columns(
        self, columns: dict[str, str]
    ) -> dict[str, str]:
        mappings: dict[str, str] = {}
        for feature_name, notion_type in columns.items():
            record_field = self._record_field_mapping_for_notion_type(
                notion_type
            ) or self._record_field_mapping_for_property_name(feature_name)
            if record_field is not None:
                mappings[feature_name] = record_field
        return mappings

    def _relation_property_matches_target_names(
        self, property_name: str, target_names: list[str]
    ) -> bool:
        property_tokens = {
            token.strip().lower().replace("_", " ")
            for token in self._strict_name_candidates(property_name)
            if isinstance(token, str) and token.strip()
        }
        target_tokens: set[str] = set()
        for target_name in target_names:
            target_tokens.update(
                token.strip().lower().replace("_", " ")
                for token in self._strict_name_candidates(target_name)
                if isinstance(token, str) and token.strip()
            )
        return len(property_tokens & target_tokens) > 0

    def _infer_notion_backward_relation_features(
        self,
        schema_spec: dict[str, dict[str, Any]],
        features_by_name: dict[str, Any],
        current_type_name: str | None = None,
        locked_feature_names: set[str] | None = None,
    ) -> dict[str, Any]:
        """Infer schema feature name -> source feature for dual Notion relation pairs."""
        locked_feature_names = locked_feature_names or set()
        logger.important(
            "backward relation: infer start "
            f"relation_props={[name for name, spec in schema_spec.items() if spec.get('type') == 'relation']} "
            f"local_features={sorted(features_by_name.keys())}"
        )
        matches: list[tuple[str, Any, str]] = []
        for property_name, property_spec in schema_spec.items():
            if property_spec.get("type") != "relation":
                continue
            if property_name in locked_feature_names:
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=already-configured"
                )
                continue
            local_feature = features_by_name.get(property_name)
            if local_feature is None:
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=local-feature-missing"
                )
                continue
            dual = property_spec.get("dual")
            if not isinstance(dual, dict):
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=dual-missing"
                )
                continue
            synced_property_name = dual.get("synced_property_name")
            target = property_spec.get("target")
            if (
                not isinstance(synced_property_name, str)
                or not synced_property_name.strip()
                or not isinstance(target, str)
                or not target.strip()
            ):
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=dual-or-target-invalid "
                    f"synced={synced_property_name!r} target={target!r}"
                )
                continue
            target_names = self._relation_target_name_candidates(target)
            is_self_referential_target = isinstance(current_type_name, str) and any(
                isinstance(target_name, str)
                and target_name.strip().lower() == current_type_name.strip().lower()
                for target_name in target_names
            )
            if (
                not self._relation_property_matches_target_names(
                    property_name, target_names
                )
                and not is_self_referential_target
            ):
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=target-name-mismatch "
                    f"target_names={target_names}"
                )
                continue
            if is_self_referential_target:
                logger.important(
                    "backward relation: self-referential target-name override "
                    f"relation={property_name!r} current_type={current_type_name!r} "
                    f"target_names={target_names}"
                )
            target_type = self._resolve_record_type_by_name_candidates(target_names)
            if target_type is None:
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=target-type-unresolved "
                    f"target_names={target_names}"
                )
                continue
            if not self._relation_feature_matches_target_type(
                local_feature, target_type
            ):
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=local-feature-target-mismatch "
                    f"local_dtype={getattr(local_feature, '_dtype_str', None)!r} "
                    f"target_type={getattr(target_type, 'name', target_type)!r}"
                )
                continue
            source_feature = None
            target_feature_type = ln.Feature.filter(
                name=target_type.name, is_type=True
            ).one_or_none()
            if target_feature_type is not None:
                source_feature = ln.Feature.filter(
                    name__iexact=synced_property_name.strip(),
                    type=target_feature_type,
                ).one_or_none()
            if (
                source_feature is None
                and getattr(target_type, "schema", None) is not None
            ):
                members = target_type.schema.members
                if hasattr(members, "filter"):
                    source_candidates = list(
                        members.filter(name__iexact=synced_property_name.strip())
                    )
                else:
                    source_candidates = [
                        feature
                        for feature in members
                        if feature.name.lower() == synced_property_name.strip().lower()
                    ]
                source_feature = self._pick_unique(source_candidates)
            if source_feature is None:
                logger.important(
                    "backward relation: skip relation "
                    f"{property_name!r} reason=source-feature-unresolved "
                    f"target_type={getattr(target_type, 'name', target_type)!r} "
                    f"synced_property={synced_property_name.strip()!r}"
                )
                continue
            logger.important(
                "backward relation: candidate "
                f"{property_name!r} -> source_feature={source_feature.name!r} "
                f"uid={source_feature.uid!r}"
            )
            matches.append(
                (property_name, source_feature, synced_property_name.strip())
            )
        if not matches:
            logger.important(
                "backward relation: infer no-candidates "
                f"matches={[(name, feature.uid) for name, feature, _ in matches]}"
            )
            return {}
        if sys.stdin is None or not sys.stdin.isatty():
            logger.warning(
                "backward relation: unresolved candidates found in "
                "a non-interactive session; skipping auto-configuration"
            )
            logger.important(
                "backward relation: infer non-interactive candidates "
                f"{[(name, feature.uid) for name, feature, _ in matches]}"
            )
            return {}

        selected_mappings: dict[str, Any] = {}
        RICH_CONSOLE.print(
            "[bold yellow]notion backward relation[/] confirm backward-derived features:",
            markup=True,
            highlight=False,
        )
        for property_name, source_feature, synced_property_name in matches:
            while True:
                selected = input(
                    f"Mark '{property_name}' as backward-derived from "
                    f"'{synced_property_name}' (source uid: {source_feature.uid})? [y/N]: "
                ).strip()
                if selected == "" or selected.lower() in {"n", "no"}:
                    logger.important(
                        "backward relation: user declined candidate "
                        f"feature={property_name!r} source_uid={source_feature.uid!r}"
                    )
                    break
                if selected.lower() in {"y", "yes"}:
                    selected_mappings[property_name] = source_feature
                    logger.important(
                        "backward relation: user accepted candidate "
                        f"feature={property_name!r} source_uid={source_feature.uid!r}"
                    )
                    break
                RICH_CONSOLE.print(
                    "[yellow]Invalid choice. Use y/yes, n/no, or Enter to skip.[/]",
                    markup=True,
                    highlight=False,
                )
        logger.important(
            "backward relation: infer confirmed-candidates "
            f"{[(name, feature.uid) for name, feature in selected_mappings.items()]}"
        )
        return selected_mappings

    def _relation_feature_matches_target_type(
        self, local_feature: Any, target_type: Any
    ) -> bool:
        """Whether a local relation feature points to a specific target record type."""
        from lamindb.models.feature import parse_dtype

        dtype_str = getattr(local_feature, "_dtype_str", None)
        if not isinstance(dtype_str, str) or not dtype_str:
            return True
        parsed = parse_dtype(dtype_str)
        if len(parsed) != 1:
            logger.important(
                "backward relation: relation-target-check "
                f"result=False reason=parsed-dtype-not-singular local_dtype={dtype_str!r}"
            )
            return False
        parsed_dtype = parsed[0]
        if parsed_dtype.get("registry_str") != "Record":
            logger.important(
                "backward relation: relation-target-check "
                f"result=False reason=registry-not-record local_dtype={dtype_str!r} "
                f"registry_str={parsed_dtype.get('registry_str')!r}"
            )
            return False
        parsed_type_uid = parsed_dtype.get("type_uid")
        target_uid = getattr(target_type, "uid", None)
        if isinstance(parsed_type_uid, str) and isinstance(target_uid, str):
            matched = parsed_type_uid == target_uid
            logger.important(
                "backward relation: relation-target-check "
                f"result={matched} reason=type-uid "
                f"parsed_type_uid={parsed_type_uid!r} target_uid={target_uid!r}"
            )
            return matched
        registry = parsed_dtype.get("registry")
        registry_uid = getattr(registry, "uid", None)
        registry_name = getattr(registry, "name", None)
        target_name = getattr(target_type, "name", None)
        if not isinstance(registry_name, str) and isinstance(registry_uid, str):
            registry_record = ln.Record.filter(
                uid=registry_uid, is_type=True
            ).one_or_none()
            registry_name = getattr(registry_record, "name", None)
        if not isinstance(target_name, str) and isinstance(target_uid, str):
            target_record = ln.Record.filter(uid=target_uid, is_type=True).one_or_none()
            target_name = getattr(target_record, "name", None)
        if isinstance(registry_name, str) and isinstance(target_name, str):
            matched = registry_name.lower() == target_name.lower()
            logger.important(
                "backward relation: relation-target-check "
                f"result={matched} reason=exact-name "
                f"registry_name={registry_name!r} target_name={target_name!r} "
                f"registry_uid={registry_uid!r} target_uid={target_uid!r}"
            )
            return matched
        if isinstance(registry_uid, str) and isinstance(target_uid, str):
            matched = registry_uid == target_uid
            logger.important(
                "backward relation: relation-target-check "
                f"result={matched} reason=uid-fallback "
                f"registry_uid={registry_uid!r} target_uid={target_uid!r} "
                f"registry_name={registry_name!r} target_name={target_name!r}"
            )
            return matched
        logger.important(
            "backward relation: relation-target-check "
            f"result=False reason=insufficient-identity "
            f"registry_uid={registry_uid!r} target_uid={target_uid!r} "
            f"registry_name={registry_name!r} target_name={target_name!r}"
        )
        return False

    @staticmethod
    def _index_feature_name_from_columns(columns: dict[str, str]) -> str | None:
        for name, notion_type in columns.items():
            if notion_type == "title":
                return name
        return None

    def _database_feature_plan(
        self,
        database_id: str,
        db_name: str | None = None,
        columns: dict[str, str] | None = None,
        schema_spec: dict[str, dict[str, Any]] | None = None,
        *,
        apply: bool = False,
        report: SyncReport | None = None,
    ) -> list[tuple[str, str, Any]]:
        db_name = db_name or database_id
        if schema_spec is None:
            schema_spec = self.reader.schema(database_id)
        if columns is None:
            columns = {k: v["type"] for k, v in schema_spec.items()}
        ordered_feature_names = list(columns)
        plan: list[tuple[str, str, Any]] = []
        for name in ordered_feature_names:
            property_spec = schema_spec.get(name, {"type": columns[name]})
            dtype_label, dtype = self._dtype_from_notion_property(
                db_name,
                name,
                property_spec,
                apply=apply,
                report=report,
            )
            plan.append(
                (
                    name,
                    dtype_label,
                    dtype,
                )
            )
        return plan

    @staticmethod
    def _append_unique(values: list[str], value: str) -> None:
        if value not in values:
            values.append(value)

    @staticmethod
    def _feature_plan_detail(
        db_name: str,
        feature_name: str,
        dtype_label: str,
        dtype: Any,
        record_field_mappings: dict[str, str],
        *,
        index_feature_name: str | None = None,
    ) -> str:
        detail = f"{db_name} / {feature_name}: {dtype_label}"
        mapped_targets: list[str] = []
        if index_feature_name is not None and feature_name == index_feature_name:
            mapped_targets.append("Schema.index")
        mapped_field = record_field_mappings.get(feature_name)
        if mapped_field is not None:
            mapped_targets.append(f"Record.{mapped_field}")
        if mapped_targets:
            detail = f"{detail} -> {' / '.join(mapped_targets)}"
        registry_targets: list[str] = []
        dtype_args = getattr(dtype, "__args__", ())
        if getattr(dtype, "__origin__", None) is list and dtype_args:
            registry_type = dtype_args[0]
            if registry_type is ln.Project:
                registry_targets.append("LaminDB.Project registry")
            elif registry_type is ln.Reference:
                registry_targets.append("LaminDB.Reference registry")
        if registry_targets:
            suffix = " / ".join(registry_targets)
            detail = (
                f"{detail} -> {suffix}"
                if " -> " not in detail
                else f"{detail} / {suffix}"
            )
        return detail

    @staticmethod
    def _serialize_feature_dtype(dtype: Any) -> str | None:
        from lamindb.models.feature import serialize_dtype

        if isinstance(dtype, str):
            return dtype
        try:
            return serialize_dtype(dtype)
        except Exception:
            return None

    def _can_safely_upgrade_feature_dtype(
        self, feature: Any, current_dtype_str: str, planned_dtype_str: str
    ) -> bool:
        from lamindb.models.feature import parse_dtype

        # Notion people properties used to be created as list[str]; upgrade to
        # list[User] so re-syncing an existing schema can recover.
        if current_dtype_str == "list[str]" and planned_dtype_str in {
            "list[cat[User]]",
            "list[User]",
        }:
            return True
        try:
            current = parse_dtype(current_dtype_str)
            planned = parse_dtype(planned_dtype_str)
        except Exception:
            return False
        if len(current) != 1 or len(planned) != 1:
            return False
        current_component = current[0]
        planned_component = planned[0]
        current_registry = current_component.get("registry_str")
        planned_registry = planned_component.get("registry_str")
        if bool(current_component.get("list")) != bool(planned_component.get("list")):
            return False
        if current_component.get("filter_str") not in {"", None}:
            return False
        if current_registry == planned_registry:
            return (
                current_component.get("type_uid") is None
                and planned_component.get("type_uid") is not None
            )
        if (
            current_registry == "Record"
            and planned_registry in {"Reference", "Project"}
            and isinstance(current_component.get("type_uid"), str)
            and planned_component.get("type_uid") in {None, ""}
            and planned_component.get("filter_str") in {"", None}
        ):
            if ln.models.RecordRecord.filter(feature=feature).exists():
                return False
            source_type = ln.Record.filter(
                uid=current_component["type_uid"], is_type=True
            ).one_or_none()
            if source_type is None:
                # Stale typed-Record dtype with no stored values: safe to retarget.
                return True
            source_name = (source_type.name or "").strip().lower()
            target_name = planned_registry.lower()
            if source_name not in {target_name, f"{target_name}s"}:
                return False
            return True
        return False

    @staticmethod
    def _record_type_move_detail(rec_type: Any, parent_type: Any) -> str:
        previous_parent_raw = getattr(getattr(rec_type, "type", None), "name", None)
        previous_parent = (
            previous_parent_raw
            if isinstance(previous_parent_raw, str) and previous_parent_raw
            else None
        )
        previous = previous_parent if previous_parent else "<root>"
        target = getattr(parent_type, "name", None) or str(parent_type)
        return f"{rec_type.name}: {previous} -> {target}"

    def _plan_or_create_db_metadata(
        self,
        db_name: str,
        feature_plan: list[tuple[str, str, Any]],
        schema_spec: dict[str, dict[str, Any]] | None = None,
        index_feature_name: str | None = None,
        record_field_mappings: dict[str, str] | None = None,
        *,
        apply: bool,
        report: SyncReport,
    ) -> tuple[Any, list[Any], Any]:
        record_field_mappings = dict(record_field_mappings or {})
        if (
            index_feature_name is not None
            and record_field_mappings.get(index_feature_name) == "name"
        ):
            # The index feature is already persisted on Record.name automatically.
            record_field_mappings.pop(index_feature_name, None)
        feature_type_qs = ln.Feature.filter(name=db_name, is_type=True)
        feature_type_count = feature_type_qs.count()
        if feature_type_count > 1:
            raise ValueError(
                f"Ambiguous LaminDB feature type name {db_name!r}: found {feature_type_count} matches."
            )
        feature_type = feature_type_qs.one_or_none()
        schema_qs = ln.Schema.filter(name=db_name)
        schema_count = schema_qs.count()
        if schema_count > 1:
            raise ValueError(
                f"Ambiguous LaminDB schema name {db_name!r}: found {schema_count} matches."
            )
        schema = schema_qs.one_or_none()
        logger.important(
            f"notion sync metadata: db={db_name!r}, apply={apply}, "
            f"feature_plan_size={len(feature_plan)}, feature_type_exists={feature_type is not None}"
        )
        if feature_type is None:
            if apply:
                feature_type = ln.Feature(name=db_name, is_type=True).save()
                self._append_unique(report.created_feature_types, db_name)
                logger.important(
                    f"notion sync metadata: created feature type {db_name!r}"
                )
            else:
                self._append_unique(report.create_feature_types, db_name)

        feature_names = [name for name, _, _ in feature_plan]
        feature_plan_by_name = {
            name: (dtype_label, dtype) for name, dtype_label, dtype in feature_plan
        }
        schema_members: list[Any] = []
        schema_members_by_name: dict[str, Any] = {}
        if schema is not None:
            members = schema.members
            schema_members = (
                list(members.all()) if hasattr(members, "all") else list(members)
            )
            schema_members_by_name = {
                feature.name: feature for feature in schema_members
            }

        existing_names = set(schema_members_by_name)
        missing_specs: list[tuple[str, str, Any]] = [
            spec for spec in feature_plan if spec[0] not in existing_names
        ]
        logger.important(
            f"notion sync metadata: db={db_name!r}, existing_features={len(existing_names)}, "
            f"missing_features={len(missing_specs)}"
        )

        type_update_specs: list[tuple[str, str, Any, Any]] = []
        dtype_update_specs: list[tuple[str, str, Any, Any, str]] = []
        if schema is not None:
            for name, dtype_label, dtype in feature_plan:
                existing_feature = schema_members_by_name.get(name)
                if existing_feature is None:
                    continue
                existing_dtype = getattr(
                    existing_feature, "_dtype_str", None
                ) or getattr(existing_feature, "dtype_as_str", None)
                planned_dtype_str = self._serialize_feature_dtype(dtype)
                if feature_type is None:
                    type_update_specs.append(
                        (name, dtype_label, dtype, existing_feature)
                    )
                    continue
                if getattr(existing_feature, "type_id", None) != getattr(
                    feature_type, "id", None
                ):
                    type_update_specs.append(
                        (name, dtype_label, dtype, existing_feature)
                    )
                if (
                    isinstance(existing_dtype, str)
                    and isinstance(planned_dtype_str, str)
                    and existing_dtype != planned_dtype_str
                ):
                    dtype_update_specs.append(
                        (
                            name,
                            dtype_label,
                            dtype,
                            existing_feature,
                            planned_dtype_str,
                        )
                    )

        if missing_specs:
            for name, dtype_label, dtype in missing_specs:
                detail = self._feature_plan_detail(
                    db_name,
                    name,
                    dtype_label,
                    dtype,
                    record_field_mappings,
                    index_feature_name=index_feature_name,
                )
                if apply:
                    self._append_unique(report.created_features, detail)
                else:
                    self._append_unique(report.create_features, detail)

        feature_update_specs: dict[str, tuple[str, str, Any]] = {}
        for name, dtype_label, dtype, _ in type_update_specs:
            feature_update_specs[name] = (name, dtype_label, dtype)
        for name, dtype_label, dtype, _, _ in dtype_update_specs:
            feature_update_specs[name] = (name, dtype_label, dtype)

        if not apply and feature_update_specs:
            for name, dtype_label, dtype in feature_update_specs.values():
                detail = self._feature_plan_detail(
                    db_name,
                    name,
                    dtype_label,
                    dtype,
                    record_field_mappings,
                    index_feature_name=index_feature_name,
                )
                self._append_unique(report.update_features, detail)

        if apply and feature_type is not None and missing_specs:
            for name, _, dtype in missing_specs:
                mapped_field = record_field_mappings.get(name)
                feature_kwargs: dict[str, Any] = {
                    "name": name,
                    "dtype": dtype,
                    "type": feature_type,
                }
                if mapped_field is not None:
                    feature_kwargs["values_through"] = mapped_field
                ln.Feature(**feature_kwargs).save()
            logger.important(
                f"notion sync metadata: created {len(missing_specs)} features for {db_name!r}"
            )

        type_updates_applied_features: set[str] = set()
        if (
            apply
            and feature_type is not None
            and schema is not None
            and type_update_specs
        ):
            for name, _, _, _ in type_update_specs:
                feature = schema_members_by_name.get(name)
                if feature is None:
                    continue
                if getattr(feature, "type_id", None) != getattr(
                    feature_type, "id", None
                ):
                    feature.type = feature_type
                    feature.save(update_fields=["type"])
                    type_updates_applied_features.add(name)
            logger.important(
                f"notion sync metadata: updated {len(type_update_specs)} features to type {db_name!r}"
            )

        dtype_updates_applied_features: set[str] = set()
        if apply and schema is not None and dtype_update_specs:
            dtype_updates_applied = 0
            for name, _, _, feature, planned_dtype_str in dtype_update_specs:
                current_dtype = getattr(feature, "_dtype_str", None)
                if not isinstance(current_dtype, str):
                    continue
                if current_dtype == planned_dtype_str:
                    continue
                if not self._can_safely_upgrade_feature_dtype(
                    feature, current_dtype, planned_dtype_str
                ):
                    logger.warning(
                        "notion sync metadata: skipped dtype update "
                        f"for feature {name!r} in {db_name!r}; "
                        f"current={current_dtype!r}, planned={planned_dtype_str!r}"
                    )
                    continue
                feature._dtype_str = planned_dtype_str
                feature.save(update_fields=["_dtype_str"])
                dtype_updates_applied += 1
                dtype_updates_applied_features.add(name)
            if dtype_updates_applied > 0:
                logger.important(
                    f"notion sync metadata: updated {dtype_updates_applied} feature dtypes for {db_name!r}"
                )
        if apply:
            applied_updates = (
                type_updates_applied_features | dtype_updates_applied_features
            )
            for name, dtype_label, dtype in feature_update_specs.values():
                if name not in applied_updates:
                    continue
                detail = self._feature_plan_detail(
                    db_name,
                    name,
                    dtype_label,
                    dtype,
                    record_field_mappings,
                    index_feature_name=index_feature_name,
                )
                self._append_unique(report.updated_features, detail)

        if feature_type is not None:
            features = list(
                ln.Feature.filter(name__in=feature_names, type=feature_type)
            )
        elif schema is not None:
            members = schema.members
            if hasattr(members, "filter"):
                features = list(members.filter(name__in=feature_names))
            else:
                features = [
                    feature for feature in members if feature.name in feature_names
                ]
        else:
            features = []
        features_by_name = {feature.name: feature for feature in features}
        locked_feature_names: set[str] = set()
        if schema is not None:
            for feature_name, feature in features_by_name.items():
                source_uid = None
                values_through_feature = getattr(feature, "values_through", None)
                values_through_uid = getattr(values_through_feature, "uid", None)
                if isinstance(values_through_uid, str):
                    source_uid = values_through_uid
                else:
                    aux = getattr(feature, "_aux", None)
                    if isinstance(aux, dict):
                        aux_source_uid = aux.get("vf")
                        if isinstance(aux_source_uid, str):
                            source_uid = aux_source_uid
                if source_uid is not None:
                    locked_feature_names.add(feature_name)

        backward_features_by_name: dict[str, Any] = {}
        if schema_spec:
            backward_features_by_name = self._infer_notion_backward_relation_features(
                schema_spec,
                features_by_name,
                current_type_name=db_name,
                locked_feature_names=locked_feature_names,
            )
        backward_feature_updates: dict[str, Any] = {}
        for target_feature_name, source_feature in backward_features_by_name.items():
            target_feature = features_by_name.get(target_feature_name)
            source_uid = getattr(source_feature, "uid", None)
            if target_feature is None or not isinstance(source_uid, str):
                continue
            existing_source_uid = None
            existing_values_through = getattr(target_feature, "values_through", None)
            existing_values_through_uid = getattr(existing_values_through, "uid", None)
            if isinstance(existing_values_through_uid, str):
                existing_source_uid = existing_values_through_uid
            else:
                target_aux = getattr(target_feature, "_aux", None)
                if isinstance(target_aux, dict):
                    aux_source_uid = target_aux.get("vf")
                    if isinstance(aux_source_uid, str):
                        existing_source_uid = aux_source_uid
            if existing_source_uid != source_uid:
                backward_feature_updates[target_feature_name] = source_feature

        if backward_feature_updates:
            for target_feature_name in sorted(backward_feature_updates):
                target_feature = features_by_name[target_feature_name]
                dtype_label, dtype = feature_plan_by_name.get(
                    target_feature.name,
                    (str(getattr(target_feature, "dtype_as_str", "unknown")), str),
                )
                detail = self._feature_plan_detail(
                    db_name,
                    target_feature.name,
                    dtype_label,
                    dtype,
                    record_field_mappings,
                    index_feature_name=index_feature_name,
                )
                if apply:
                    self._append_unique(report.updated_features, detail)
                else:
                    self._append_unique(report.update_features, detail)
        if apply:
            for target_feature_name, source_feature in backward_feature_updates.items():
                target_feature = features_by_name[target_feature_name]
                target_feature.values_through = source_feature
                target_feature.save()

        if schema is None:
            if apply:
                index_feature = next(
                    (
                        feature
                        for feature in features
                        if feature.name == index_feature_name
                    ),
                    None,
                )
                logger.important(
                    f"notion sync metadata: creating schema {db_name!r} with "
                    f"{len(features)} features, index={index_feature_name!r}"
                )
                schema_features: list[Any] = []
                for feature in features:
                    if index_feature is not None and feature is index_feature:
                        # schema.index persists on Record.name automatically and the
                        # index feature is attached through the dedicated index field.
                        continue
                    mapped_field = record_field_mappings.get(feature.name)
                    if (
                        mapped_field is not None
                        and getattr(feature, "values_through", None) != mapped_field
                    ):
                        feature.values_through = mapped_field
                        feature.save()
                    schema_features.append(feature)
                schema = ln.Schema(
                    schema_features,
                    name=db_name,
                    index=index_feature,
                    # Notion properties are all optional; a sparse set_values()
                    # would fail COLUMN_NOT_IN_DATAFRAME if members were required.
                    minimal_set=False,
                ).save()
                self._append_unique(report.created_schemas, db_name)
                logger.important(f"notion sync metadata: created schema {db_name!r}")
            else:
                self._append_unique(report.create_schemas, db_name)
        else:
            logger.important(f"notion sync metadata: schema {db_name!r} already exists")
            if getattr(schema, "minimal_set", True) is True:
                if apply:
                    schema.minimal_set = False
                    schema.save()
                    self._append_unique(report.updated_schemas, db_name)
                    logger.important(
                        f"notion sync metadata: set minimal_set=False on schema {db_name!r}"
                    )
                else:
                    self._append_unique(report.update_schemas, db_name)
            if missing_specs:
                if apply:
                    self._append_unique(report.updated_schemas, db_name)
                else:
                    self._append_unique(report.update_schemas, db_name)
            if apply and record_field_mappings:
                for feature in features:
                    mapped_field = record_field_mappings.get(feature.name)
                    if mapped_field is None:
                        continue
                    if getattr(feature, "values_through", None) != mapped_field:
                        feature.values_through = mapped_field
                        feature.save()
                        logger.important(
                            "notion sync metadata: set values_through="
                            f"{mapped_field!r} on feature {feature.name!r}"
                        )
        return feature_type, features, schema

    def _create_record_type(
        self,
        database_id: str,
        db_name: str,
        db_description: str | None,
        db_emoji: str | None,
        report: SyncReport,
        parent_type=None,
    ):
        schema_spec = self.reader.schema(database_id)
        columns = {k: v["type"] for k, v in schema_spec.items()}
        feature_plan = self._database_feature_plan(
            database_id,
            db_name=db_name,
            columns=columns,
            schema_spec=schema_spec,
            apply=True,
            report=report,
        )
        index_feature_name = self._index_feature_name_from_columns(columns)
        record_field_mappings = self._record_field_mappings_from_columns(columns)
        _, _, schema = self._plan_or_create_db_metadata(
            db_name,
            feature_plan,
            schema_spec=schema_spec,
            index_feature_name=index_feature_name,
            record_field_mappings=record_field_mappings,
            apply=True,
            report=report,
        )
        assert schema is not None  # schema is always created/resolved in apply mode
        record_kwargs: dict[str, Any] = {
            "name": db_name,
            "description": db_description,
            "is_type": True,
            "schema": schema,
        }
        if parent_type is not None:
            record_kwargs["type"] = parent_type
        aux = self._merge_aux_with_emoji(None, db_emoji)
        if aux is not None:
            record_kwargs["_aux"] = aux
        return ln.Record(**record_kwargs).save()

    def _resolve_record_type(
        self,
        database_id: str,
        *,
        apply: bool,
        report: SyncReport,
        plan_schema: bool = True,
        payload: dict | None = None,
        parent_type=None,
        parent_types_by_page_id: dict[str, Any] | None = None,
    ):
        if payload is None:
            payload = self.reader._call("GET", f"/databases/{database_id}")
        if parent_type is None and parent_types_by_page_id:
            parent_page_id = self._database_parent_pages.get(
                _normalize_notion_id(database_id) or database_id
            )
            if parent_page_id is not None:
                parent_type = parent_types_by_page_id.get(parent_page_id)
        if parent_type is None and parent_types_by_page_id:
            parent_page_id = self._database_parent_page_id(payload)
            if parent_page_id is not None:
                parent_type = parent_types_by_page_id.get(parent_page_id)
        if parent_type is None and len(parent_types_by_page_id or {}) == 1:
            parent_type = next(iter(parent_types_by_page_id.values()))
        db_name = self._database_title(payload, fallback=database_id)
        db_description = self._database_description(payload)
        db_emoji = self._database_emoji(payload)
        if db_emoji is None and parent_type is not None:
            parent_aux = getattr(parent_type, "_aux", None)
            if isinstance(parent_aux, dict):
                db_emoji = parent_aux.get("ei")
        qs = ln.Record.filter(name=db_name, is_type=True)
        count = qs.count()
        if not plan_schema:
            if count == 0:
                if not apply:
                    report.create_record_types.append(db_name)
                    return None
                rec_type_kwargs: dict[str, Any] = {
                    "name": db_name,
                    "description": db_description,
                    "is_type": True,
                }
                if parent_type is not None:
                    rec_type_kwargs["type"] = parent_type
                aux = self._merge_aux_with_emoji(None, db_emoji)
                if aux is not None:
                    rec_type_kwargs["_aux"] = aux
                rec_type = ln.Record(**rec_type_kwargs).save()
                report.created_record_types.append(db_name)
                return rec_type
            if count > 1:
                raise ValueError(
                    f"Ambiguous Lamin record type name {db_name!r}: found {count} matches."
                )
            rec_type = qs.one()
            if apply:
                changed = False
                if rec_type.description != db_description:
                    rec_type.description = db_description
                    changed = True
                merged_aux = self._merge_aux_with_emoji(
                    getattr(rec_type, "_aux", None), db_emoji
                )
                if getattr(rec_type, "_aux", None) != merged_aux:
                    rec_type._aux = merged_aux
                    changed = True
                if parent_type is not None and getattr(
                    rec_type, "type_id", None
                ) != getattr(parent_type, "id", None):
                    rec_type.type = parent_type
                    changed = True
                    detail = self._record_type_move_detail(rec_type, parent_type)
                    self._append_unique(report.updated_record_types, detail)
                if changed:
                    rec_type.save()
            elif parent_type is not None and getattr(
                rec_type, "type_id", None
            ) != getattr(parent_type, "id", None):
                detail = self._record_type_move_detail(rec_type, parent_type)
                self._append_unique(report.update_record_types, detail)
            return rec_type
        if count == 0:
            schema_spec = self.reader.schema(database_id)
            relation_debug = {
                name: {
                    "target": spec.get("target"),
                    "synced_property_name": (
                        spec.get("dual", {}).get("synced_property_name")
                        if isinstance(spec.get("dual"), dict)
                        else None
                    ),
                }
                for name, spec in schema_spec.items()
                if spec.get("type") == "relation"
            }
            logger.important(
                "relation sync: loaded schema spec "
                f"db={db_name!r} relation_props={relation_debug}"
            )
            columns = {k: v["type"] for k, v in schema_spec.items()}
            feature_plan = self._database_feature_plan(
                database_id,
                db_name=db_name,
                columns=columns,
                schema_spec=schema_spec,
                apply=apply,
                report=report,
            )
            index_feature_name = self._index_feature_name_from_columns(columns)
            record_field_mappings = self._record_field_mappings_from_columns(columns)
            if not apply:
                report.create_record_types.append(db_name)
                self._plan_or_create_db_metadata(
                    db_name,
                    feature_plan,
                    schema_spec=schema_spec,
                    index_feature_name=index_feature_name,
                    record_field_mappings=record_field_mappings,
                    apply=False,
                    report=report,
                )
                return None
            if parent_type is not None:
                rec_type = self._create_record_type(
                    database_id,
                    db_name,
                    db_description,
                    db_emoji,
                    report=report,
                    parent_type=parent_type,
                )
            else:
                rec_type = self._create_record_type(
                    database_id, db_name, db_description, db_emoji, report=report
                )
            report.created_record_types.append(db_name)
            return rec_type
        if count > 1:
            raise ValueError(
                f"Ambiguous Lamin record type name {db_name!r}: found {count} matches."
            )
        rec_type = qs.one()
        schema_spec = self.reader.schema(database_id)
        relation_debug = {
            name: {
                "target": spec.get("target"),
                "synced_property_name": (
                    spec.get("dual", {}).get("synced_property_name")
                    if isinstance(spec.get("dual"), dict)
                    else None
                ),
            }
            for name, spec in schema_spec.items()
            if spec.get("type") == "relation"
        }
        logger.important(
            "relation sync: loaded schema spec "
            f"db={db_name!r} relation_props={relation_debug}"
        )
        columns = {k: v["type"] for k, v in schema_spec.items()}
        feature_plan = self._database_feature_plan(
            database_id,
            db_name=db_name,
            columns=columns,
            schema_spec=schema_spec,
            apply=apply,
            report=report,
        )
        index_feature_name = self._index_feature_name_from_columns(columns)
        record_field_mappings = self._record_field_mappings_from_columns(columns)
        record_type_move_detail: str | None = None
        if parent_type is not None and getattr(rec_type, "type_id", None) != getattr(
            parent_type, "id", None
        ):
            record_type_move_detail = self._record_type_move_detail(
                rec_type, parent_type
            )
        if apply:
            _, features, schema = self._plan_or_create_db_metadata(
                db_name,
                feature_plan,
                schema_spec=schema_spec,
                index_feature_name=index_feature_name,
                record_field_mappings=record_field_mappings,
                apply=True,
                report=report,
            )
            changed = False
            if rec_type.description != db_description:
                rec_type.description = db_description
                changed = True
            merged_aux = self._merge_aux_with_emoji(
                getattr(rec_type, "_aux", None), db_emoji
            )
            if getattr(rec_type, "_aux", None) != merged_aux:
                rec_type._aux = merged_aux
                changed = True
            if rec_type.schema is None and schema is not None:
                rec_type.schema = schema
                changed = True
            if parent_type is not None and getattr(
                rec_type, "type_id", None
            ) != getattr(parent_type, "id", None):
                rec_type.type = parent_type
                changed = True
                if record_type_move_detail is not None:
                    self._append_unique(
                        report.updated_record_types, record_type_move_detail
                    )
            if rec_type.schema is not None:
                existing_schema_feature_names = {
                    feature.name for feature in rec_type.schema.members
                }
                missing_schema_features = [
                    feature
                    for feature in features
                    if feature.name not in existing_schema_feature_names
                ]
                if missing_schema_features:
                    rec_type.schema.add(missing_schema_features)
            if changed:
                rec_type.save()
        else:
            self._plan_or_create_db_metadata(
                db_name,
                feature_plan,
                schema_spec=schema_spec,
                index_feature_name=index_feature_name,
                record_field_mappings=record_field_mappings,
                apply=False,
                report=report,
            )
            if rec_type.schema is None:
                self._append_unique(
                    report.update_record_types,
                    self._record_type_schema_attach_detail(rec_type.name, db_name),
                )
            if record_type_move_detail is not None:
                self._append_unique(report.update_record_types, record_type_move_detail)
        return rec_type

    def _validate_schema(
        self, database_id: str, rec_type, *, apply: bool = True
    ) -> None:
        notion_props = set(self.reader.columns(database_id))
        if rec_type.schema is None and not apply:
            logger.important(
                f"notion sync schema-check: discovered schemaless record type "
                f"{rec_type.name!r}; dry run reports planned schema attachment"
            )
            return
        schema_features = self._schema_feature_names(rec_type)
        missing_features = sorted(notion_props - schema_features)
        extra_features = sorted(schema_features - notion_props)
        if missing_features and not apply and not extra_features:
            logger.important(
                f"notion sync schema-check: discovered existing schema {rec_type.name!r} "
                f"with missing features={missing_features}; dry run reports planned updates"
            )
            return
        if missing_features or extra_features:
            problems: list[str] = []
            if missing_features:
                problems.append(f"missing in Lamin schema: {missing_features}")
            if extra_features:
                problems.append(f"extra in Lamin schema: {extra_features}")
            msg = "; ".join(problems)
            raise ValueError(
                f"Schema mismatch for database {_compact_uuid(database_id)!r} <-> record type "
                f"{rec_type.name!r}: {msg}"
            )

    @staticmethod
    def _existing_edit_map(by_id: dict) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for notion_id, record in by_id.items():
            out[notion_id] = _normalized_timestamp(getattr(record, "updated_at", None))
        return out

    def import_pages(
        self,
        parents: str | list[str],
        *,
        apply: bool = False,
        depth: int | None = None,
    ) -> SyncReport:
        """Import parent trees, validating schema before any write.

        `parents` are Notion page/database IDs. The sync discovers databases under
        these roots, validates property parity against Lamin schemas, then performs
        an idempotent upsert/materialize pass.
        """
        if isinstance(parents, str):
            parent_ids = [parents]
        else:
            parent_ids = list(parents)
        if not parent_ids:
            raise ValueError("parents is required and must contain at least one ID.")

        report = SyncReport(
            apply=apply,
            message="Dry run report -- nothing got created" if not apply else None,
        )
        db_id_set, parent_pages = self._collect_database_ids(parent_ids, depth=depth)
        db_ids = sorted(db_id_set)

        # Parent pages can also map to LaminDB record types.
        parent_types_by_page_id: dict[str, Any] = {}
        for parent_id, parent_type_name in sorted(parent_pages.items()):
            parent_type = self._resolve_or_create_type_by_name(
                parent_type_name,
                emoji=self._parent_page_emojis.get(parent_id),
                apply=apply,
                report=report,
            )
            if parent_type is not None:
                parent_types_by_page_id[parent_id] = parent_type
        for child_page_id, ancestor_page_id in sorted(
            self._parent_page_parents.items()
        ):
            child_type = parent_types_by_page_id.get(child_page_id)
            ancestor_type = parent_types_by_page_id.get(ancestor_page_id)
            if child_type is None or ancestor_type is None:
                continue
            if getattr(child_type, "type_id", None) == getattr(
                ancestor_type, "id", None
            ):
                continue
            detail = self._record_type_move_detail(child_type, ancestor_type)
            if apply:
                child_type.type = ancestor_type
                child_type.save(update_fields=["type"])
                self._append_unique(report.updated_record_types, detail)
            else:
                self._append_unique(report.update_record_types, detail)

        if not db_ids:
            report.discovered_pages = len(parent_pages)
            if depth == 0:
                return report
            raise ValueError(
                "No child databases discovered under parents. In phase 1, sync operates "
                "on page trees that include at least one Notion database."
            )
        report.databases = [_compact_uuid(db_id) for db_id in db_ids]

        # Step 1: resolve and validate schema parity before any write.
        rec_types: dict[str, Any] = {}
        db_specs: dict[str, dict[str, Any]] = {}
        for db_id in db_ids:
            logger.important(
                f"notion sync schema-check: resolving record type for db={_compact_uuid(db_id)}"
            )
            normalized_db_id = _normalize_notion_id(db_id) or db_id
            seed_page_ids = self._seed_page_ids_by_database.get(normalized_db_id, set())
            rec_type = self._resolve_record_type(
                db_id,
                apply=apply,
                report=report,
                plan_schema=not bool(seed_page_ids),
                parent_types_by_page_id=parent_types_by_page_id,
            )
            if rec_type is not None:
                self._validate_schema(db_id, rec_type, apply=apply)
                logger.important(
                    f"notion sync schema-check: validated db={_compact_uuid(db_id)} against "
                    f"record_type={rec_type.name!r}"
                )
            rec_types[db_id] = rec_type
            db_specs[db_id] = self.reader.schema(db_id)

        after_maps: dict[str, dict[str, Any]] = {}
        to_write: dict[str, list[dict[str, Any]]] = {}
        planned_transfers_by_db: dict[str, dict[str, str]] = {}

        with _bulk_creation():
            # Phase A: discover + upsert identity rows.
            for db_id in db_ids:
                rec_type = rec_types[db_id]
                normalized_db_id = _normalize_notion_id(db_id) or db_id
                seed_page_ids = self._seed_page_ids_by_database.get(
                    normalized_db_id, set()
                )
                if seed_page_ids:
                    rows = self._rows_for_seed_pages(
                        db_id, seed_page_ids, include_page_emoji=apply
                    )
                else:
                    rows = self.reader.rows(db_id, include_page_emoji=apply)
                report.discovered += len(rows)
                logger.important(
                    f"notion sync phase A: db={_compact_uuid(db_id)}, discovered_rows={len(rows)}"
                )
                if rec_type is None:
                    # dry-run mode with missing type: all discovered rows are new.
                    report.created += len(rows)
                    to_write[db_id] = []
                    after_maps[db_id] = {}
                    _, _, file_props = _kinds(db_specs[db_id])
                    transfer_map, transfer_details = _planned_missing_file_transfers(
                        rows, file_props
                    )
                    note_transfer_map, note_transfer_details = (
                        _planned_missing_embedded_transfers(self.reader, rows)
                        if not apply
                        else ({}, [])
                    )
                    transfer_map.update(note_transfer_map)
                    transfer_details.extend(note_transfer_details)
                    planned_transfers_by_db[db_id] = transfer_map
                    for detail in transfer_details:
                        if detail not in report.create_artifacts:
                            report.create_artifacts.append(detail)
                    continue

                before = _existing_by_ref(rec_type)
                before_edit = self._existing_edit_map(before)

                writes: list[dict[str, Any]] = []
                force_materialize_seed_rows = apply and bool(seed_page_ids)
                preview_seed_rows_in_dry_run = (not apply) and bool(seed_page_ids)
                for row in rows:
                    notion_id = row["notion_id"]
                    edited = _normalized_timestamp(
                        _parse_notion_timestamp(row.get("last_edited_time"))
                    )
                    existing = before_edit.get(notion_id)
                    if notion_id not in before:
                        report.created += 1
                        writes.append(row)
                    elif existing == edited:
                        report.unchanged += 1
                        if force_materialize_seed_rows:
                            # Seed-page sync needs idempotent feature/notes backfill
                            # even when timestamps are unchanged (e.g. prior partial runs).
                            writes.append(row)
                    else:
                        report.updated += 1
                        writes.append(row)
                to_write[db_id] = writes
                preview_rows = rows if preview_seed_rows_in_dry_run else writes
                logger.important(
                    f"notion sync phase A: db={_compact_uuid(db_id)}, create={sum(1 for row in writes if row['notion_id'] not in before)}, "
                    f"update={sum(1 for row in writes if row['notion_id'] in before)}, unchanged={len(rows) - len(writes)}"
                )
                _, _, file_props = _kinds(db_specs[db_id])
                transfer_map, transfer_details = _planned_missing_file_transfers(
                    preview_rows, file_props
                )
                if not apply:
                    note_transfer_map, note_transfer_details = (
                        _planned_missing_embedded_transfers(self.reader, preview_rows)
                    )
                    transfer_map.update(note_transfer_map)
                    transfer_details.extend(note_transfer_details)
                planned_transfers_by_db[db_id] = transfer_map
                if not apply:
                    for detail in transfer_details:
                        if detail not in report.create_artifacts:
                            report.create_artifacts.append(detail)
                    rel, _, _ = _kinds(db_specs[db_id])
                    if rel and preview_rows:
                        feat = _feat_map(rec_type.schema)
                        _, relation_pending = _resolve_relation_records_for_rows(
                            self.reader,
                            preview_rows,
                            rel,
                            feat,
                            None,
                            rec_type=rec_type,
                            apply=False,
                            report=report,
                        )
                        report.pending_relations += relation_pending

                if apply:
                    logger.important(
                        f"notion sync phase A: upserting identity rows for db={_compact_uuid(db_id)}"
                    )
                    after_maps[db_id] = _upsert_all(rec_type, rows)
                else:
                    after_maps[db_id] = before

            report.discovered_pages = (
                len(parent_pages) + len(report.databases) + report.discovered
            )
            if not apply:
                return report

            # Phase B: materialize only changed/new rows.
            for db_id in db_ids:
                rec_type = rec_types[db_id]
                spec = db_specs[db_id]
                write_rows = to_write[db_id]
                if not write_rows:
                    continue
                logger.important(
                    f"notion sync phase B: materializing db={_compact_uuid(db_id)} rows={len(write_rows)}"
                )
                stats = _write(
                    self.reader,
                    write_rows,
                    rec_type,
                    spec,
                    prop_map=None,
                    by_id=after_maps[db_id],
                    transfer_details_by_url=planned_transfers_by_db.get(db_id, {}),
                    report=report,
                )
                report.pending_relations += stats["pending"]
                logger.important(
                    f"notion sync phase B: finished db={_compact_uuid(db_id)} records={stats['records']} pending_relations={stats['pending']}"
                )

        logger.important(
            f"notion sync done: created={report.created}, updated={report.updated}, unchanged={report.unchanged}"
        )
        return report


class RecordSyncer(_NotionSyncer):
    """Record-specific Notion sync utility."""


class ProjectSyncer:
    """Project-specific Notion sync utility."""

    SKIP_PROJECT_RECORD_PROPERTIES: set[str] = {"presentations", "meetings"}

    PROJECT_STATUS_TO_CODE: dict[str, int] = {
        str(status): code for status, code in LAMIN_PROJECT_STATUS_TO_CODE.items()
    }
    PROJECT_STATUS_ALIASES: dict[str, str] = {
        "on hold": "paused",
        "on_hold": "paused",
        "in progress": "active",
        "in_progress": "active",
        "up next": "up-next",
        "up_next": "up-next",
        "complete": "completed",
        "completed": "completed",
        "continued": "background",
        "cancelled": "canceled",
    }

    def __init__(self, reader: _NotionReader) -> None:
        self.reader = reader

    @staticmethod
    def _normalize_name(value: str) -> str:
        return value.strip().lower().replace("-", " ").replace("_", " ")

    def _normalize_status(self, value: str) -> str:
        normalized = self._normalize_name(value)
        return self.PROJECT_STATUS_ALIASES.get(normalized, normalized)

    def _status_mapping_error(
        self, db_name: str, unknown_statuses: list[str], known_statuses: list[str]
    ) -> ValueError:
        unknown = ", ".join(sorted(set(unknown_statuses)))
        expected = ", ".join(sorted(set(known_statuses)))
        return ValueError(
            "Project status mapping mismatch for Notion database "
            f"{db_name!r}. Unknown Notion status labels: {unknown}. "
            "Please update status names in Notion to match LaminDB status names: "
            f"{expected}."
        )

    @staticmethod
    def _is_project_name(name: str) -> bool:
        normalized = name.strip().lower()
        return normalized in {"project", "projects", "task", "tasks"}

    @staticmethod
    def _is_reference_name(name: str) -> bool:
        normalized = name.strip().lower()
        return normalized in {"reference", "references"}

    @staticmethod
    def _is_user_directory_name(name: str) -> bool:
        normalized = name.strip().lower()
        return normalized in {
            "user",
            "users",
            "person",
            "people",
            "members",
            "team members",
        }

    @staticmethod
    def _is_parent_relation_name(name: str) -> bool:
        normalized = name.strip().lower().replace("_", " ")
        parent_tokens = {
            "parent",
            "parents",
            "project",
            "program",
            "initiative",
            "portfolio",
            "owner project",
            "superproject",
        }
        return any(token in normalized for token in parent_tokens)

    @staticmethod
    def _is_child_relation_name(name: str) -> bool:
        normalized = name.strip().lower().replace("_", " ")
        child_tokens = {
            "child",
            "children",
            "task",
            "tasks",
            "subproject",
            "sub project",
            "sub-project",
        }
        return any(token in normalized for token in child_tokens)

    @staticmethod
    def _is_predecessor_relation_name(name: str) -> bool:
        normalized = name.strip().lower().replace("_", " ")
        predecessor_tokens = {
            "predecessor",
            "predecessors",
            "dependency",
            "dependencies",
            "blocked by",
            "depends on",
            "requires",
        }
        return any(token in normalized for token in predecessor_tokens)

    @staticmethod
    def _is_successor_relation_name(name: str) -> bool:
        normalized = name.strip().lower().replace("_", " ")
        successor_tokens = {
            "successor",
            "successors",
            "dependent",
            "dependents",
            "follow up",
            "follow-up",
            "follows",
            "after",
        }
        return any(token in normalized for token in successor_tokens)

    @staticmethod
    def _notion_project_url(notion_id: str) -> str:
        compact_id = _normalize_notion_id(notion_id) or notion_id
        return f"https://notion.so/laminlabs/{compact_id}"

    @staticmethod
    def _append_unmapped(
        *,
        report: SyncReport,
        db_name: str,
        property_name: str,
        notion_type: str,
        reason: str,
    ) -> None:
        detail = f"{db_name} / {property_name} ({notion_type}): {reason}"
        _append_unique(report.unmapped_properties, detail)
        logger.warning(f"notion project sync unmapped property: {detail}")

    def build_mapping(
        self,
        *,
        db_name: str,
        schema_spec: dict[str, dict[str, Any]],
        report: SyncReport,
    ) -> dict[str, Any]:
        mapping: dict[str, Any] = {
            "title": None,
            "description": None,
            "timeline": None,
            "start_date": None,
            "end_date": None,
            "status": None,
            "created_time": None,
            "last_edited_time": None,
            "created_by": None,
            "parents_rel": set(),
            "children_rel": set(),
            "predecessors_rel": set(),
            "successors_rel": set(),
            "references_rel": set(),
            "record_rel": {},
            "people_roles": {},
        }
        for property_name, property_spec in schema_spec.items():
            notion_type = property_spec.get("type")
            normalized = self._normalize_name(property_name)
            if notion_type == "title":
                if mapping["title"] is None:
                    mapping["title"] = property_name
                else:
                    self._append_unmapped(
                        report=report,
                        db_name=db_name,
                        property_name=property_name,
                        notion_type=notion_type,
                        reason="extra title field (only one title can map to Project.name)",
                    )
                continue
            if notion_type == "rich_text" and normalized in {"summary", "description"}:
                if mapping["description"] is None:
                    mapping["description"] = property_name
                else:
                    self._append_unmapped(
                        report=report,
                        db_name=db_name,
                        property_name=property_name,
                        notion_type=notion_type,
                        reason="duplicate description-like field",
                    )
                continue
            if notion_type == "date":
                if normalized == "timeline":
                    mapping["timeline"] = property_name
                    continue
                if normalized in {"start", "start date", "start_date"}:
                    mapping["start_date"] = property_name
                    continue
                if normalized in {
                    "end",
                    "end date",
                    "end_date",
                    "deadline",
                    "due",
                    "due date",
                    "due_date",
                }:
                    mapping["end_date"] = property_name
                    continue
                self._append_unmapped(
                    report=report,
                    db_name=db_name,
                    property_name=property_name,
                    notion_type=notion_type,
                    reason="ambiguous date field name for Project start/end mapping",
                )
                continue
            if notion_type in {"status", "select"} and normalized == "status":
                mapping["status"] = property_name
                continue
            if notion_type == "created_time":
                mapping["created_time"] = property_name
                continue
            if notion_type == "last_edited_time":
                mapping["last_edited_time"] = property_name
                continue
            if notion_type == "created_by":
                mapping["created_by"] = property_name
                continue
            if notion_type == "people":
                role = normalized.replace("_", " ").strip() or "member"
                mapping["people_roles"][property_name] = role
                continue
            if notion_type == "relation":
                target = property_spec.get("target")
                target_names = []
                if isinstance(target, str) and target:
                    target_names = self._target_names(target)
                if normalized == "responsible":
                    mapping["people_roles"][property_name] = "responsible"
                    continue
                if any(self._is_project_name(name) for name in target_names):
                    if self._is_predecessor_relation_name(property_name):
                        mapping["predecessors_rel"].add(property_name)
                    elif self._is_successor_relation_name(property_name):
                        mapping["successors_rel"].add(property_name)
                    elif self._is_parent_relation_name(property_name):
                        mapping["parents_rel"].add(property_name)
                    elif self._is_child_relation_name(property_name):
                        mapping["children_rel"].add(property_name)
                    else:
                        self._append_unmapped(
                            report=report,
                            db_name=db_name,
                            property_name=property_name,
                            notion_type=notion_type,
                            reason=(
                                "project relation has no recognized semantic "
                                "(expected parent/child or predecessor/successor semantics)"
                            ),
                        )
                    continue
                if any(self._is_reference_name(name) for name in target_names):
                    mapping["references_rel"].add(property_name)
                    continue
                record_type = self._resolve_record_type_by_name_candidates(target_names)
                if record_type is not None:
                    normalized_property_name = (
                        property_name.strip().lower().replace(" ", "_")
                    )
                    if normalized_property_name in self.SKIP_PROJECT_RECORD_PROPERTIES:
                        self._append_unmapped(
                            report=report,
                            db_name=db_name,
                            property_name=property_name,
                            notion_type=notion_type,
                            reason=(
                                "intentionally skipped ProjectRecord mapping "
                                "(populate from Record side)"
                            ),
                        )
                        continue
                    mapping["record_rel"][property_name] = record_type
                    target_name = getattr(record_type, "name", None)
                    if not isinstance(target_name, str) or not target_name:
                        target_name = getattr(record_type.__class__, "__name__", None)
                    if not isinstance(target_name, str) or not target_name:
                        target_name = str(record_type)
                    _append_unique(
                        report.mapped_project_record_relations,
                        f"{db_name} / {property_name} -> ProjectRecord(feature={property_name}, target={target_name})",
                    )
                    continue
                self._append_unmapped(
                    report=report,
                    db_name=db_name,
                    property_name=property_name,
                    notion_type=notion_type,
                    reason="relation target does not map to Project/Reference/Record type",
                )
                continue
            self._append_unmapped(
                report=report,
                db_name=db_name,
                property_name=property_name,
                notion_type=str(notion_type),
                reason="unsupported Project field mapping",
            )
        return mapping

    def _target_names(self, target: str) -> list[str]:
        db_payload = self._safe_call(f"/databases/{target}")
        if db_payload is not None:
            db_name = _NotionSyncer._database_title(db_payload, fallback=target)
            return [db_name]
        data_source_payload = self._safe_call(f"/data_sources/{target}")
        if data_source_payload is None:
            return []
        names: list[str] = []
        data_source_name = data_source_payload.get("name")
        if isinstance(data_source_name, str) and data_source_name.strip():
            names.append(data_source_name.strip())
        parent = data_source_payload.get("parent")
        if isinstance(parent, dict):
            parent_db_id = parent.get("database_id")
            if isinstance(parent_db_id, str) and parent_db_id.strip():
                parent_db_payload = self._safe_call(f"/databases/{parent_db_id}")
                if parent_db_payload is not None:
                    names.append(
                        _NotionSyncer._database_title(
                            parent_db_payload, fallback=parent_db_id
                        )
                    )
        return list(dict.fromkeys(names))

    def _safe_call(self, path: str) -> dict | None:
        try:
            return self.reader._call("GET", path)
        except LookupError:
            return None
        except httpx.HTTPStatusError as error:
            if error.response is not None and error.response.status_code == 400:
                return None
            raise

    @staticmethod
    def _pick_unique(records: list[Any]) -> Any | None:
        unique_by_id: dict[Any, Any] = {}
        for record in records:
            record_id = getattr(record, "id", None)
            key = record_id if record_id is not None else id(record)
            unique_by_id[key] = record
        if len(unique_by_id) == 1:
            return next(iter(unique_by_id.values()))
        return None

    def _resolve_record_type_by_name_candidates(self, names: list[str]) -> Any | None:
        matches: list[Any] = []
        for candidate in names:
            if not isinstance(candidate, str) or not candidate.strip():
                continue
            qs = ln.Record.filter(name__iexact=candidate.strip(), is_type=True)
            if qs.count() == 1:
                matches.append(qs.one())
        return self._pick_unique(matches)

    def validate_status_mapping(
        self,
        *,
        db_name: str,
        schema_spec: dict[str, dict[str, Any]],
        mapping: dict[str, Any],
        rows: list[dict[str, Any]] | None = None,
    ) -> None:
        status_property = mapping.get("status")
        if status_property is None:
            return
        known = list(self.PROJECT_STATUS_TO_CODE.keys())
        status_spec = schema_spec.get(status_property, {})
        choices = status_spec.get("choices")
        unknown_from_schema: list[str] = []
        if isinstance(choices, list):
            for choice in choices:
                if not isinstance(choice, str) or not choice.strip():
                    continue
                if self._normalize_status(choice) not in self.PROJECT_STATUS_TO_CODE:
                    unknown_from_schema.append(choice)
        if unknown_from_schema:
            raise self._status_mapping_error(db_name, unknown_from_schema, known)
        if rows:
            unknown_from_rows: list[str] = []
            for row in rows:
                value = row.get(status_property)
                if not isinstance(value, str) or not value.strip():
                    continue
                if self._normalize_status(value) not in self.PROJECT_STATUS_TO_CODE:
                    unknown_from_rows.append(value)
            if unknown_from_rows:
                raise self._status_mapping_error(db_name, unknown_from_rows, known)

    @staticmethod
    def _existing_edit_map(by_id: dict[str, Any]) -> dict[str, datetime | None]:
        out: dict[str, datetime | None] = {}
        for notion_id, project in by_id.items():
            out[notion_id] = _normalized_timestamp(getattr(project, "updated_at", None))
        return out

    def existing_by_notion_id(self, rows: list[dict[str, Any]]) -> dict[str, Any]:
        notion_ids = {
            _normalize_notion_id(row.get("notion_id"))
            for row in rows
            if isinstance(row.get("notion_id"), str) and row.get("notion_id")
        }
        notion_ids = {notion_id for notion_id in notion_ids if notion_id is not None}
        if not notion_ids:
            return {}
        url_to_id = {
            self._notion_project_url(notion_id): notion_id for notion_id in notion_ids
        }
        by_id: dict[str, Any] = {}
        for project in ln.Project.filter(url__in=list(url_to_id)):
            notion_id = url_to_id.get(getattr(project, "url", None))
            if notion_id is not None and notion_id in notion_ids:
                by_id[notion_id] = project
        return by_id

    def upsert_all(
        self,
        *,
        rows: list[dict[str, Any]],
        by_id: dict[str, Any],
        title_property: str | None,
        project_type: Any | None = None,
    ) -> dict[str, Any]:
        for row in rows:
            raw_notion_id = row.get("notion_id")
            if not isinstance(raw_notion_id, str) or not raw_notion_id:
                continue
            notion_id = _normalize_notion_id(raw_notion_id) or raw_notion_id
            project_url = self._notion_project_url(notion_id)
            row_emoji = row.get("__notion_emoji__")
            created_at = _parse_notion_timestamp(row.get("created_time"))
            updated_at = _parse_notion_timestamp(row.get("last_edited_time"))
            if created_at is None:
                created_at = updated_at
            if updated_at is None:
                updated_at = created_at
            row_name = row.get(title_property) if title_property else row.get("name")
            project = by_id.get(notion_id)
            if project is None:
                project_kwargs: dict[str, Any] = {
                    "name": row_name or _compact_uuid(notion_id),
                    "url": project_url,
                }
                if project_type is not None:
                    project_kwargs["type"] = project_type
                project = ln.Project(**project_kwargs)
                project._aux = _NotionSyncer._merge_aux_with_emoji(None, row_emoji)
                if created_at is not None:
                    project.created_at = created_at
                if updated_at is not None:
                    project.updated_at = updated_at
                by_id[notion_id] = project.save()
                continue
            changed_fields: list[str] = []
            if getattr(project, "url", None) != project_url:
                project.url = project_url
                changed_fields.append("url")
            if isinstance(row_name, str) and row_name and row_name != project.name:
                project.name = row_name
                changed_fields.append("name")
            if project_type is not None and getattr(
                project, "type_id", None
            ) != getattr(project_type, "id", None):
                project.type = project_type
                changed_fields.append("type")
            if created_at is not None and _normalized_timestamp(
                project.created_at
            ) != _normalized_timestamp(created_at):
                project.created_at = created_at
                changed_fields.append("created_at")
            if updated_at is not None and _normalized_timestamp(
                project.updated_at
            ) != _normalized_timestamp(updated_at):
                project.updated_at = updated_at
                changed_fields.append("updated_at")
            merged_aux = _NotionSyncer._merge_aux_with_emoji(
                getattr(project, "_aux", None), row_emoji
            )
            if getattr(project, "_aux", None) != merged_aux:
                project._aux = merged_aux
                changed_fields.append("_aux")
            if changed_fields:
                project.save(update_fields=changed_fields)
        return by_id

    def _date_from_iso(self, value: Any) -> date | None:
        if not isinstance(value, str) or not value:
            return None
        parsed = _parse_notion_timestamp(value)
        if parsed is not None:
            return parsed.date()
        try:
            return date.fromisoformat(value)
        except ValueError:
            return None

    def _date_range_from_page(
        self, page_payload: dict[str, Any], property_name: str
    ) -> tuple[date | None, date | None]:
        properties = page_payload.get("properties")
        if not isinstance(properties, dict):
            return None, None
        property_payload = properties.get(property_name)
        if not isinstance(property_payload, dict):
            return None, None
        if property_payload.get("type") != "date":
            return None, None
        date_payload = property_payload.get("date")
        if not isinstance(date_payload, dict):
            return None, None
        start_date = self._date_from_iso(date_payload.get("start"))
        end_date = self._date_from_iso(date_payload.get("end"))
        return start_date, end_date

    def _status_code_from_row(
        self, db_name: str, row: dict[str, Any], status_property: str | None
    ) -> int | None:
        if status_property is None:
            return None
        raw_value = row.get(status_property)
        if not isinstance(raw_value, str) or not raw_value.strip():
            return None
        normalized = self._normalize_status(raw_value)
        if normalized not in self.PROJECT_STATUS_TO_CODE:
            raise self._status_mapping_error(
                db_name, [raw_value], list(self.PROJECT_STATUS_TO_CODE)
            )
        return self.PROJECT_STATUS_TO_CODE[normalized]

    def _sync_project_user_role(
        self, project: Any, role: str, users: list[Any]
    ) -> None:
        desired_ids = {getattr(user, "id", None) for user in users}
        desired_ids.discard(None)
        existing_links = list(project.links_user.filter(role=role))
        existing_by_user_id = {
            getattr(link, "user_id", None): link
            for link in existing_links
            if getattr(link, "user_id", None) is not None
        }
        for user_id, link in existing_by_user_id.items():
            if user_id not in desired_ids:
                link.delete()
        existing_ids = set(existing_by_user_id)
        for user in users:
            user_id = getattr(user, "id", None)
            if user_id is None or user_id in existing_ids:
                continue
            project.links_user.create(user=user, role=role)

    def write_projects(
        self,
        *,
        db_name: str,
        rows: list[dict[str, Any]],
        by_id: dict[str, Any],
        mapping: dict[str, Any],
        report: SyncReport | None = None,
        transfer_details_by_url: dict[str, str] | None = None,
    ) -> dict[str, int]:
        artifacts_by_url = _batch_artifacts(
            rows,
            set(),
            transfer_details_by_url=transfer_details_by_url,
            report=report,
        )
        markdown_by_notion_id: dict[str, str] = {}
        embedded_transfer_details_by_url: dict[str, str] = {}
        embedded_file_urls: set[str] = set()
        for row in rows:
            notion_id = row.get("notion_id")
            if not isinstance(notion_id, str) or not notion_id:
                continue
            markdown_content = self.reader.page_markdown(notion_id)
            markdown_by_notion_id[notion_id] = markdown_content
            per_page_transfer_map, _ = _planned_embedded_file_transfers(
                markdown_content, notion_id
            )
            for url, detail in per_page_transfer_map.items():
                embedded_file_urls.add(url)
                embedded_transfer_details_by_url.setdefault(url, detail)
        missing_embedded_urls = embedded_file_urls - set(artifacts_by_url)
        if missing_embedded_urls:
            artifacts_by_url.update(
                _ensure_artifacts(
                    missing_embedded_urls,
                    transfer_details_by_url=embedded_transfer_details_by_url,
                    report=report,
                    with_key=False,
                    kind="__easset__",
                    description="imported from Notion",
                )
            )

        relation_props = (
            set(mapping["parents_rel"])
            | set(mapping["children_rel"])
            | set(mapping["predecessors_rel"])
            | set(mapping["successors_rel"])
            | set(mapping["references_rel"])
            | set(mapping["record_rel"].keys())
        )
        fake_features: dict[str, Any] = {}
        for prop in (
            mapping["parents_rel"]
            | mapping["children_rel"]
            | mapping["predecessors_rel"]
            | mapping["successors_rel"]
        ):
            fake_features[prop] = type(
                "FeatureProxy", (), {"name": prop, "_dtype_str": "list[cat[Project]]"}
            )()
        for prop in mapping["references_rel"]:
            fake_features[prop] = type(
                "FeatureProxy", (), {"name": prop, "_dtype_str": "list[cat[Reference]]"}
            )()
        for prop, target_record_type in mapping["record_rel"].items():
            target_uid = getattr(target_record_type, "uid", None)
            if not isinstance(target_uid, str) or not target_uid:
                continue
            fake_features[prop] = type(
                "FeatureProxy",
                (),
                {"name": prop, "_dtype_str": f"list[cat[Record[{target_uid}]]]"},
            )()
        if relation_props:
            resolved_relations, relation_pending = _resolve_relation_records_for_rows(
                self.reader,
                rows,
                relation_props,
                fake_features,
                None,
                rec_type=type("ProjectRegistry", (), {"name": db_name})(),
                apply=True,
                report=report,
            )
        else:
            resolved_relations, relation_pending = {}, 0

        people_link_stats: dict[str, dict[str, set[str]]] = {
            prop: {"resolved": set(), "unresolved": set()}
            for prop in mapping["people_roles"]
        }
        user_relation_ids: set[str] = set()
        for people_property in mapping["people_roles"]:
            for row in rows:
                value = row.get(people_property)
                values = value if isinstance(value, list) else [value]
                for item in values:
                    if isinstance(item, str) and item:
                        user_relation_ids.add(_normalize_notion_id(item) or item)
        created_by_property = mapping.get("created_by")
        for row in rows:
            created_by_value = (
                row.get(created_by_property) if created_by_property else None
            )
            if isinstance(created_by_value, str) and created_by_value:
                user_relation_ids.add(
                    _normalize_notion_id(created_by_value) or created_by_value
                )
        users_by_notion_id = _resolved_users_by_notion_id(
            self.reader, sorted(user_relation_ids)
        )
        record_link_features: dict[str, Any] = {}
        if mapping["record_rel"]:
            feature_type = ln.Feature.filter(name=db_name, is_type=True).one_or_none()
            if feature_type is None:
                feature_type = ln.Feature(name=db_name, is_type=True).save()
            for prop_name, target_record_type in mapping["record_rel"].items():
                relation_dtype = _NotionSyncer._list_dtype_for(target_record_type)
                feature = ln.Feature.filter(
                    name__iexact=prop_name,
                    type=feature_type,
                ).one_or_none()
                if feature is None:
                    feature = ln.Feature(
                        name=prop_name,
                        dtype=relation_dtype,
                        type=feature_type,
                    ).save()
                record_link_features[prop_name] = feature

        records = 0
        for row in rows:
            raw_notion_id = row.get("notion_id")
            if not isinstance(raw_notion_id, str) or not raw_notion_id:
                continue
            notion_id = _normalize_notion_id(raw_notion_id) or raw_notion_id
            project = by_id.get(notion_id)
            if project is None:
                continue
            changed_fields: list[str] = []
            title_property = mapping.get("title")
            if title_property:
                title = row.get(title_property)
                if isinstance(title, str) and title and title != project.name:
                    project.name = title
                    changed_fields.append("name")
            description_property = mapping.get("description")
            if description_property:
                description = row.get(description_property)
                if description != getattr(project, "description", None):
                    project.description = description
                    changed_fields.append("description")
            status_code = self._status_code_from_row(
                db_name, row, mapping.get("status")
            )
            if status_code is not None and status_code != getattr(
                project, "_status_code", None
            ):
                project._status_code = status_code
                changed_fields.append("_status_code")
            created_by_property = mapping.get("created_by")
            if created_by_property:
                created_by_notion = row.get(created_by_property)
                if isinstance(created_by_notion, str):
                    normalized_user_id = (
                        _normalize_notion_id(created_by_notion) or created_by_notion
                    )
                    resolved_user = users_by_notion_id.get(normalized_user_id)
                    if resolved_user is not None and getattr(
                        project, "created_by_id", None
                    ) != getattr(resolved_user, "id", None):
                        project.created_by = resolved_user
                        changed_fields.append("created_by")
            created_at = _parse_notion_timestamp(row.get("created_time"))
            updated_at = _parse_notion_timestamp(row.get("last_edited_time"))
            if created_at is not None and _normalized_timestamp(
                getattr(project, "created_at", None)
            ) != _normalized_timestamp(created_at):
                project.created_at = created_at
                changed_fields.append("created_at")
            if updated_at is not None and _normalized_timestamp(
                getattr(project, "updated_at", None)
            ) != _normalized_timestamp(updated_at):
                project.updated_at = updated_at
                changed_fields.append("updated_at")

            start_date: date | None = None
            end_date: date | None = None
            start_property = mapping.get("start_date")
            if start_property:
                start_date = self._date_from_iso(row.get(start_property))
            end_property = mapping.get("end_date")
            if end_property:
                end_date = self._date_from_iso(row.get(end_property))
            timeline_property = mapping.get("timeline")
            if timeline_property and (start_date is None or end_date is None):
                page_payload = self.reader._call(
                    "GET", f"/pages/{_notion_api_id(notion_id)}"
                )
                timeline_start, timeline_end = self._date_range_from_page(
                    page_payload, timeline_property
                )
                if start_date is None:
                    start_date = timeline_start
                if end_date is None:
                    end_date = timeline_end
            if start_date != getattr(project, "start_date", None):
                project.start_date = start_date
                changed_fields.append("start_date")
            if end_date != getattr(project, "end_date", None):
                project.end_date = end_date
                changed_fields.append("end_date")

            if changed_fields:
                project.save(update_fields=changed_fields)

            def _resolved_values(
                row_values: dict[str, Any], prop_names: set[str], model_name: str
            ) -> list[Any]:
                notion_ids: list[str] = []
                for prop_name in prop_names:
                    value = row_values.get(prop_name)
                    raw_values = value if isinstance(value, list) else [value]
                    for item in raw_values:
                        if isinstance(item, str) and item:
                            notion_ids.append(_normalize_notion_id(item) or item)
                values = [
                    resolved_relations[relation_id]
                    for relation_id in notion_ids
                    if relation_id in resolved_relations
                ]
                return [
                    value
                    for value in values
                    if value.__class__.__name__.lower() == model_name.lower()
                ]

            project.parents.set(
                _resolved_values(row, mapping["parents_rel"], "Project")
            )
            project.children.set(
                _resolved_values(row, mapping["children_rel"], "Project")
            )
            project.predecessors.set(
                _resolved_values(row, mapping["predecessors_rel"], "Project")
            )
            project.successors.set(
                _resolved_values(row, mapping["successors_rel"], "Project")
            )
            project.references.set(
                _resolved_values(row, mapping["references_rel"], "Reference")
            )
            for prop_name, feature in record_link_features.items():
                linked_records = _resolved_values(row, {prop_name}, "Record")
                desired_ids = {getattr(record, "id", None) for record in linked_records}
                desired_ids.discard(None)
                existing_links = list(project.links_record.filter(feature=feature))
                existing_by_record_id = {
                    getattr(link, "record_id", None): link
                    for link in existing_links
                    if getattr(link, "record_id", None) is not None
                }
                for record_id, link in existing_by_record_id.items():
                    if record_id not in desired_ids:
                        link.delete()
                existing_ids = set(existing_by_record_id)
                for record in linked_records:
                    record_id = getattr(record, "id", None)
                    if record_id is None or record_id in existing_ids:
                        continue
                    project.links_record.create(record=record, feature=feature)
            for people_property, role in mapping["people_roles"].items():
                user_ids = row.get(people_property)
                raw_ids = user_ids if isinstance(user_ids, list) else [user_ids]
                users: list[Any] = []
                stats = people_link_stats[people_property]
                for user_id in raw_ids:
                    if not isinstance(user_id, str) or not user_id:
                        continue
                    normalized_user_id = _normalize_notion_id(user_id) or user_id
                    user = users_by_notion_id.get(normalized_user_id)
                    if user is not None:
                        stats["resolved"].add(normalized_user_id)
                        users.append(user)
                    else:
                        stats["unresolved"].add(normalized_user_id)
                self._sync_project_user_role(project, role=role, users=users)

            markdown_content = markdown_by_notion_id.get(notion_id)
            if markdown_content is None:
                markdown_content = self.reader.page_markdown(notion_id)
            _attach_project_markdown(
                project,
                _rewrite_embedded_file_refs(markdown_content, artifacts_by_url),
            )
            records += 1
        if report is not None:
            for people_property, stats in sorted(people_link_stats.items()):
                resolved_count = len(stats["resolved"])
                pending_count = len(stats["unresolved"])
                relation_pending += pending_count
                unresolved_names: list[str] = []
                unresolved_lookup: list[str] = []
                for notion_id in sorted(stats["unresolved"]):
                    name, lookup_source = _notion_user_or_page_lookup(
                        self.reader, notion_id
                    )
                    display_name = name or (_compact_uuid(notion_id) or notion_id)
                    unresolved_names.append(display_name)
                    unresolved_lookup.append(f"{display_name} [{lookup_source}]")
                pending_name_note = (
                    f", pending_names={unresolved_names}" if unresolved_names else ""
                )
                pending_lookup_note = (
                    f", pending_lookup={unresolved_lookup}" if unresolved_lookup else ""
                )
                _append_unique(
                    report.relation_value_links,
                    f"{db_name} / {people_property}: "
                    f"resolved_existing={resolved_count}, "
                    "stub_create=0, "
                    f"pending_unresolved={pending_count}{pending_name_note}"
                    f"{pending_lookup_note}",
                )
        return {"records": records, "pending": relation_pending}


class NotionSyncer(RecordSyncer):
    """Dispatcher that routes Notion databases to record or project syncers."""

    def __init__(self, token: str | None = None) -> None:
        super().__init__(token=token)
        self.project_syncer = ProjectSyncer(self.reader)

    def _is_project_database(self, payload: dict[str, Any], database_id: str) -> bool:
        db_name = self._database_title(payload, fallback=database_id)
        normalized_name = db_name.strip().lower()
        if normalized_name in {"task", "tasks"}:
            return True
        resolved_type, _ = self._resolve_or_plan_relation_record_type(
            [db_name], apply=False, report=None, prefer_plural=False
        )
        return resolved_type is ln.Project

    @staticmethod
    def _project_type_name_for_database(db_name: str) -> str | None:
        normalized_name = db_name.strip().lower()
        if normalized_name in {"task", "tasks"}:
            return "Tasks"
        if normalized_name in {"project", "projects"}:
            return "Work packages"
        return None

    @staticmethod
    def _resolve_or_create_project_type(type_name: str, *, apply: bool) -> Any | None:
        qs = ln.Project.filter(name=type_name, is_type=True)
        count = qs.count()
        if count > 1:
            raise ValueError(
                f"Ambiguous LaminDB project type name {type_name!r}: found {count} matches."
            )
        if count == 1:
            return qs.one()
        if not apply:
            return None
        return ln.Project(name=type_name, is_type=True).save()

    def import_pages(
        self,
        parents: str | list[str],
        *,
        apply: bool = False,
        depth: int | None = None,
    ) -> SyncReport:
        if isinstance(parents, str):
            parent_ids = [parents]
        else:
            parent_ids = list(parents)
        if not parent_ids:
            raise ValueError("parents is required and must contain at least one ID.")

        report = SyncReport(
            apply=apply,
            message="Dry run report -- nothing got created" if not apply else None,
        )
        db_id_set, parent_pages = self._collect_database_ids(parent_ids, depth=depth)
        db_ids = sorted(db_id_set)
        if not db_ids:
            report.discovered_pages = len(parent_pages)
            if depth == 0:
                return report
            raise ValueError(
                "No child databases discovered under parents. In phase 1, sync operates "
                "on page trees that include at least one Notion database."
            )
        report.databases = [_compact_uuid(db_id) for db_id in db_ids]

        db_payloads: dict[str, dict[str, Any]] = {}
        project_db_ids: set[str] = set()
        record_db_ids: set[str] = set()
        project_type_by_db_id: dict[str, Any | None] = {}
        for db_id in db_ids:
            payload = self.reader._call("GET", f"/databases/{db_id}")
            db_payloads[db_id] = payload
            if self._is_project_database(payload, db_id):
                project_db_ids.add(db_id)
                db_name = self._database_title(payload, fallback=db_id)
                project_type_name = self._project_type_name_for_database(db_name)
                if project_type_name is not None:
                    project_type_by_db_id[db_id] = self._resolve_or_create_project_type(
                        project_type_name, apply=apply
                    )
                else:
                    project_type_by_db_id[db_id] = None
            else:
                record_db_ids.add(db_id)

        parent_types_by_page_id: dict[str, Any] = {}
        if record_db_ids:
            for parent_id, parent_type_name in sorted(parent_pages.items()):
                parent_type = self._resolve_or_create_type_by_name(
                    parent_type_name,
                    emoji=self._parent_page_emojis.get(parent_id),
                    apply=apply,
                    report=report,
                )
                if parent_type is not None:
                    parent_types_by_page_id[parent_id] = parent_type
            for child_page_id, ancestor_page_id in sorted(
                self._parent_page_parents.items()
            ):
                child_type = parent_types_by_page_id.get(child_page_id)
                ancestor_type = parent_types_by_page_id.get(ancestor_page_id)
                if child_type is None or ancestor_type is None:
                    continue
                if getattr(child_type, "type_id", None) == getattr(
                    ancestor_type, "id", None
                ):
                    continue
                detail = self._record_type_move_detail(child_type, ancestor_type)
                if apply:
                    child_type.type = ancestor_type
                    child_type.save(update_fields=["type"])
                    self._append_unique(report.updated_record_types, detail)
                else:
                    self._append_unique(report.update_record_types, detail)

        rec_types: dict[str, Any] = {}
        db_specs: dict[str, dict[str, Any]] = {}
        project_mappings: dict[str, dict[str, Any]] = {}
        for db_id in db_ids:
            schema_spec = self.reader.schema(db_id)
            db_specs[db_id] = schema_spec
            payload = db_payloads[db_id]
            db_name = self._database_title(payload, fallback=db_id)
            if db_id in project_db_ids:
                mapping = self.project_syncer.build_mapping(
                    db_name=db_name, schema_spec=schema_spec, report=report
                )
                project_mappings[db_id] = mapping
                self.project_syncer.validate_status_mapping(
                    db_name=db_name, schema_spec=schema_spec, mapping=mapping
                )
                rec_types[db_id] = None
                continue
            logger.important(
                f"notion sync schema-check: resolving record type for db={_compact_uuid(db_id)}"
            )
            normalized_db_id = _normalize_notion_id(db_id) or db_id
            seed_page_ids = self._seed_page_ids_by_database.get(normalized_db_id, set())
            rec_type = self._resolve_record_type(
                db_id,
                apply=apply,
                report=report,
                plan_schema=not bool(seed_page_ids),
                parent_types_by_page_id=parent_types_by_page_id,
            )
            if rec_type is not None:
                self._validate_schema(db_id, rec_type, apply=apply)
                logger.important(
                    f"notion sync schema-check: validated db={_compact_uuid(db_id)} against "
                    f"record_type={rec_type.name!r}"
                )
            rec_types[db_id] = rec_type

        after_maps: dict[str, dict[str, Any]] = {}
        to_write: dict[str, list[dict[str, Any]]] = {}
        planned_transfers_by_db: dict[str, dict[str, str]] = {}

        with _bulk_creation():
            for db_id in db_ids:
                normalized_db_id = _normalize_notion_id(db_id) or db_id
                seed_page_ids = self._seed_page_ids_by_database.get(
                    normalized_db_id, set()
                )
                if seed_page_ids:
                    rows = self._rows_for_seed_pages(
                        db_id, seed_page_ids, include_page_emoji=apply
                    )
                else:
                    rows = self.reader.rows(db_id, include_page_emoji=apply)
                report.discovered += len(rows)
                if db_id in project_db_ids:
                    db_name = self._database_title(db_payloads[db_id], fallback=db_id)
                    mapping = project_mappings[db_id]
                    self.project_syncer.validate_status_mapping(
                        db_name=db_name,
                        schema_spec=db_specs[db_id],
                        mapping=mapping,
                        rows=rows,
                    )
                    before = self.project_syncer.existing_by_notion_id(rows)
                    before_edit = self.project_syncer._existing_edit_map(before)
                    writes: list[dict[str, Any]] = []
                    force_materialize_seed_rows = apply and bool(seed_page_ids)
                    preview_seed_rows_in_dry_run = (not apply) and bool(seed_page_ids)
                    for row in rows:
                        notion_id = row["notion_id"]
                        edited = _normalized_timestamp(
                            _parse_notion_timestamp(row.get("last_edited_time"))
                        )
                        existing = before_edit.get(notion_id)
                        if notion_id not in before:
                            report.created_projects += 1
                            writes.append(row)
                        elif existing == edited:
                            report.unchanged_projects += 1
                            if force_materialize_seed_rows:
                                writes.append(row)
                        else:
                            report.updated_projects += 1
                            writes.append(row)
                    to_write[db_id] = writes
                    preview_rows = rows if preview_seed_rows_in_dry_run else writes
                    transfer_map, transfer_details = (
                        _planned_missing_embedded_transfers(self.reader, preview_rows)
                        if not apply
                        else ({}, [])
                    )
                    planned_transfers_by_db[db_id] = transfer_map
                    if not apply:
                        for detail in transfer_details:
                            if detail not in report.create_artifacts:
                                report.create_artifacts.append(detail)
                    if apply:
                        after_maps[db_id] = self.project_syncer.upsert_all(
                            rows=rows,
                            by_id=before,
                            title_property=mapping.get("title"),
                            project_type=project_type_by_db_id.get(db_id),
                        )
                    else:
                        after_maps[db_id] = before
                    continue

                rec_type = rec_types[db_id]
                if rec_type is None:
                    report.created += len(rows)
                    to_write[db_id] = []
                    after_maps[db_id] = {}
                    _, _, file_props = _kinds(db_specs[db_id])
                    transfer_map, transfer_details = _planned_missing_file_transfers(
                        rows, file_props
                    )
                    note_transfer_map, note_transfer_details = (
                        _planned_missing_embedded_transfers(self.reader, rows)
                        if not apply
                        else ({}, [])
                    )
                    transfer_map.update(note_transfer_map)
                    transfer_details.extend(note_transfer_details)
                    planned_transfers_by_db[db_id] = transfer_map
                    for detail in transfer_details:
                        if detail not in report.create_artifacts:
                            report.create_artifacts.append(detail)
                    continue

                before = _existing_by_ref(rec_type)
                before_edit = self._existing_edit_map(before)
                writes = []
                force_materialize_seed_rows = apply and bool(seed_page_ids)
                preview_seed_rows_in_dry_run = (not apply) and bool(seed_page_ids)
                for row in rows:
                    notion_id = row["notion_id"]
                    edited = _normalized_timestamp(
                        _parse_notion_timestamp(row.get("last_edited_time"))
                    )
                    existing = before_edit.get(notion_id)
                    if notion_id not in before:
                        report.created += 1
                        writes.append(row)
                    elif existing == edited:
                        report.unchanged += 1
                        if force_materialize_seed_rows:
                            writes.append(row)
                    else:
                        report.updated += 1
                        writes.append(row)
                to_write[db_id] = writes
                preview_rows = rows if preview_seed_rows_in_dry_run else writes
                _, _, file_props = _kinds(db_specs[db_id])
                transfer_map, transfer_details = _planned_missing_file_transfers(
                    preview_rows, file_props
                )
                if not apply:
                    note_transfer_map, note_transfer_details = (
                        _planned_missing_embedded_transfers(self.reader, preview_rows)
                    )
                    transfer_map.update(note_transfer_map)
                    transfer_details.extend(note_transfer_details)
                planned_transfers_by_db[db_id] = transfer_map
                if not apply:
                    for detail in transfer_details:
                        if detail not in report.create_artifacts:
                            report.create_artifacts.append(detail)
                    rel, _, _ = _kinds(db_specs[db_id])
                    if rel and preview_rows:
                        feat = _feat_map(rec_type.schema)
                        _, relation_pending = _resolve_relation_records_for_rows(
                            self.reader,
                            preview_rows,
                            rel,
                            feat,
                            None,
                            rec_type=rec_type,
                            apply=False,
                            report=report,
                        )
                        report.pending_relations += relation_pending
                if apply:
                    after_maps[db_id] = _upsert_all(rec_type, rows)
                else:
                    after_maps[db_id] = before

            report.discovered_pages = (
                len(parent_pages) + len(report.databases) + report.discovered
            )
            if not apply:
                return report

            for db_id in db_ids:
                write_rows = to_write[db_id]
                if not write_rows:
                    continue
                if db_id in project_db_ids:
                    db_name = self._database_title(db_payloads[db_id], fallback=db_id)
                    stats = self.project_syncer.write_projects(
                        db_name=db_name,
                        rows=write_rows,
                        by_id=after_maps[db_id],
                        mapping=project_mappings[db_id],
                        transfer_details_by_url=planned_transfers_by_db.get(db_id, {}),
                        report=report,
                    )
                    report.pending_relations += stats["pending"]
                    continue
                rec_type = rec_types[db_id]
                spec = db_specs[db_id]
                stats = _write(
                    self.reader,
                    write_rows,
                    rec_type,
                    spec,
                    prop_map=None,
                    by_id=after_maps[db_id],
                    transfer_details_by_url=planned_transfers_by_db.get(db_id, {}),
                    report=report,
                )
                report.pending_relations += stats["pending"]

        logger.important(
            "notion sync done: "
            f"record_create={report.created}, record_update={report.updated}, "
            f"record_unchanged={report.unchanged}, "
            f"project_create={report.created_projects}, "
            f"project_update={report.updated_projects}, "
            f"project_unchanged={report.unchanged_projects}"
        )
        return report


@ln.flow("Ofbk5ruuTiN2")
def sync_from_notion(
    *,
    parents: str | list[str],
    token: str | None = None,
    apply: bool = False,
    depth: int | None = None,
) -> SyncReport:
    """Sync Notion pages via the class-based sync API.

    Args:
        parents: Notion page or database ids.
        token: Notion API token. Defaults to the ``NOTION_TOKEN`` environment variable.
        apply: Write to LaminDB. By default this is a dry run.
        depth: How many levels of child pages and databases to walk.
            ``None`` walks the whole tree. ``0`` syncs only the given parents.
    """
    syncer = NotionSyncer(token=token)
    if isinstance(parents, str):
        parent_list = [parents]
    elif isinstance(parents, list):
        parent_list = list(parents)
    else:
        raise TypeError("parents must be a str or list[str].")
    report = syncer.import_pages(parents=parent_list, apply=apply, depth=depth)
    RICH_CONSOLE.print(report.to_pretty_text(), markup=True, highlight=False)
    return report


__all__ = [
    "sync_from_notion",
    "SyncReport",
    "NotionSyncer",
    "RecordSyncer",
    "ProjectSyncer",
]
