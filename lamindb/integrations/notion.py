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
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import httpx
from lamin_utils import logger
from rich.console import Console
from rich.markup import escape as rich_escape

import lamindb as ln

API_VERSION = "2026-03-11"
BASE = "https://api.notion.com/v1"
UUID_DASHED_PATTERN = re.compile(
    r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
)
RICH_CONSOLE = Console(force_terminal=True, no_color=False)


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
        lines.extend(
            [
                metric("create_records", self.created, action_color),
                metric("update_records", self.updated, action_color),
                metric("unchanged_records", self.unchanged, action_color),
                metric("pending_relations", self.pending_relations, action_color),
                metric("failed_records", self.failed, action_color),
            ]
        )
        if self.errors:
            lines.append("")
            lines.append("[bold red]Errors[/]")
            lines.extend(f"[red]- {error}[/]" for error in self.errors)
        return "\n".join(lines)


def _flatten(prop: dict) -> Any:
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
        return [x["id"] for x in (prop.get(t) or [])]
    if t in ("created_by", "last_edited_by"):
        u = prop.get(t)
        return u["id"] if u else None
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
                row[name] = _flatten(prop)
            # page-level fields win over any same-named user property
            row["notion_id"] = page.get("id")
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
            lines.append(f"{indent}> {text}".rstrip())
        elif block_type == "code":
            language = payload.get("language") if isinstance(payload, dict) else None
            fence = f"```{language}" if language else "```"
            lines.extend([f"{indent}{fence}", f"{indent}{text}", f"{indent}```"])
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
            payload = self._call("GET", f"/blocks/{block_id}/children", params=params)
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
            payload = self._call("GET", f"/blocks/{page_id}/children", params=params)
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
        raw = self._call("GET", f"/pages/{page_id}")
        row: dict[str, Any] = {
            "notion_id": None,
            "created_time": None,
            "last_edited_time": None,
        }
        for name, prop in raw.get("properties", {}).items():
            if name in drop or prop.get("type") in drop:
                continue
            row[name] = _flatten(prop)
        row["notion_id"] = raw.get("id")
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
            artifacts[url] = ln.Artifact(
                tmp_path, key=_artifact_key_from_url(url)
            ).save()
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
) -> dict[str, Any]:
    urls = {url for _, _, url in _iter_file_urls(rows, file_props)}
    return _ensure_artifacts(
        urls,
        transfer_details_by_url=transfer_details_by_url,
        report=report,
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
            uuids = val if isinstance(val, list) else [val]
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
            values[f] = val
    return values, pending


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

    targets = {u for row in rows for p in rel for u in (row.get(p) or [])}
    resolved = _resolved_map(targets)  # one query for all relation targets

    if by_id is None:
        by_id = _existing_by_ref(rec_type)

    records = pending = 0
    for row in rows:
        rec = by_id.get(row["notion_id"])
        if rec is None:  # not imported yet — nothing to write
            continue
        values, p = _row_values(
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
            _attach_page_markdown(rec, reader.page_markdown(notion_id))
        pending += p
        records += 1
    return {"records": records, "pending": pending}


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


def _upsert_all(rec_type, rows) -> dict:
    """Upsert every row of one database against a single existing-record map.

    Returns the {notion_uuid: record} map (existing + newly created), ready to
    hand to :func:`_write` so it never re-queries.
    """
    by_id = _existing_by_ref(rec_type)  # ONE query, not one per row
    for row in rows:
        nid, name = row["notion_id"], row.get("name")
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
                "GET", f"/blocks/{block_id}/children", params=params
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
        limit: int | None = None,
        discovered_children: list[int] | None = None,
    ) -> bool:
        for block in self._iter_block_children(block_id):
            if (
                limit is not None
                and discovered_children is not None
                and discovered_children[0] >= limit
            ):
                return True
            bid = block.get("id")
            if bid and bid in seen:
                continue
            if bid:
                seen.add(bid)
            if block.get("type") == "child_database" and bid:
                is_new = bid not in out
                out.add(bid)
                if is_new and discovered_children is not None:
                    discovered_children[0] += 1
                normalized_db_id = _normalize_notion_id(bid)
                if parent_page_id is not None and normalized_db_id is not None:
                    self._database_parent_pages[normalized_db_id] = parent_page_id
                if (
                    limit is not None
                    and discovered_children is not None
                    and discovered_children[0] >= limit
                ):
                    return True
            if block.get("has_children") and bid:
                reached_limit = self._collect_databases_from_block(
                    bid,
                    seen,
                    out,
                    parent_page_id=parent_page_id,
                    limit=limit,
                    discovered_children=discovered_children,
                )
                if reached_limit:
                    return True
        return False

    def _collect_database_ids(
        self, parents: list[str], limit: int | None = None
    ) -> tuple[set[str], dict[str, str]]:
        if limit is not None and limit < 0:
            raise ValueError("limit must be >= 0 when provided.")
        database_ids: set[str] = set()
        parent_pages: dict[str, str] = {}
        parent_page_emojis: dict[str, str | None] = {}
        self._database_parent_pages = {}
        self._parent_page_parents = {}
        seen_blocks: set[str] = set()
        discovered_children = [0]
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
                if limit == 0:
                    continue
                if limit is not None and discovered_children[0] >= limit:
                    continue
                for row in self.reader.rows(parent, limit=limit):
                    notion_id = row.get("notion_id")
                    if notion_id:
                        reached_limit = self._collect_databases_from_block(
                            notion_id,
                            seen_blocks,
                            database_ids,
                            parent_page_id=None,
                            limit=limit,
                            discovered_children=discovered_children,
                        )
                        if reached_limit:
                            break
                continue
            page_payload = self._safe_call(f"/pages/{parent}")
            if page_payload is None:
                raise LookupError(
                    f"Parent {_compact_uuid(parent)!r} is neither a readable database nor page."
                )
            parent_title = _page_title(page_payload).strip() or _compact_uuid(parent)
            parent_id = _normalize_notion_id(parent) or parent
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
            if limit == 0:
                continue
            reached_limit = self._collect_databases_from_block(
                parent,
                seen_blocks,
                database_ids,
                parent_page_id=parent_id,
                limit=limit,
                discovered_children=discovered_children,
            )
            if reached_limit:
                break
        self._parent_page_emojis = parent_page_emojis
        return database_ids, parent_pages

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
        if notion_type in {"multi_select", "people", "relation"}:
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
        if notion_type in {"multi_select", "people", "relation"}:
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
            return self._name_candidates(db_name)

        data_source_payload = self._safe_call(f"/data_sources/{target}")
        if data_source_payload is None:
            return []

        names: list[str] = []
        ds_name = data_source_payload.get("name")
        if isinstance(ds_name, str) and ds_name.strip():
            names.extend(self._name_candidates(ds_name.strip()))

        parent = data_source_payload.get("parent")
        if isinstance(parent, dict):
            parent_db_id = parent.get("database_id")
            if isinstance(parent_db_id, str) and parent_db_id.strip():
                parent_db_payload = self._safe_call(f"/databases/{parent_db_id}")
                if parent_db_payload is not None:
                    db_name = self._database_title(
                        parent_db_payload, fallback=parent_db_id
                    )
                    names.extend(self._name_candidates(db_name))

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
            for token in self._name_candidates(property_name)
            if isinstance(token, str) and token.strip()
        }
        target_tokens: set[str] = set()
        for target_name in target_names:
            target_tokens.update(
                token.strip().lower().replace("_", " ")
                for token in self._name_candidates(target_name)
                if isinstance(token, str) and token.strip()
            )
        return len(property_tokens & target_tokens) > 0

    def _infer_notion_backward_relation_feature(
        self,
        schema_spec: dict[str, dict[str, Any]],
        features_by_name: dict[str, Any],
    ) -> tuple[str | None, Any | None]:
        """Infer (feature_name, source_feature) for a dual Notion relation pair."""
        logger.important(
            "notion sync backward-debug: infer start "
            f"relation_props={[name for name, spec in schema_spec.items() if spec.get('type') == 'relation']} "
            f"local_features={sorted(features_by_name.keys())}"
        )
        matches: list[tuple[str, Any, str]] = []
        for property_name, property_spec in schema_spec.items():
            if property_spec.get("type") != "relation":
                continue
            local_feature = features_by_name.get(property_name)
            if local_feature is None:
                logger.important(
                    "notion sync backward-debug: skip relation "
                    f"{property_name!r} reason=local-feature-missing"
                )
                continue
            dual = property_spec.get("dual")
            if not isinstance(dual, dict):
                logger.important(
                    "notion sync backward-debug: skip relation "
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
                    "notion sync backward-debug: skip relation "
                    f"{property_name!r} reason=dual-or-target-invalid "
                    f"synced={synced_property_name!r} target={target!r}"
                )
                continue
            target_names = self._relation_target_name_candidates(target)
            if not self._relation_property_matches_target_names(
                property_name, target_names
            ):
                logger.important(
                    "notion sync backward-debug: skip relation "
                    f"{property_name!r} reason=target-name-mismatch "
                    f"target_names={target_names}"
                )
                continue
            target_type = self._resolve_record_type_by_name_candidates(target_names)
            if target_type is None:
                logger.important(
                    "notion sync backward-debug: skip relation "
                    f"{property_name!r} reason=target-type-unresolved "
                    f"target_names={target_names}"
                )
                continue
            if not self._relation_feature_matches_target_type(
                local_feature, target_type
            ):
                logger.important(
                    "notion sync backward-debug: skip relation "
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
                    "notion sync backward-debug: skip relation "
                    f"{property_name!r} reason=source-feature-unresolved "
                    f"target_type={getattr(target_type, 'name', target_type)!r} "
                    f"synced_property={synced_property_name.strip()!r}"
                )
                continue
            logger.important(
                "notion sync backward-debug: candidate relation "
                f"{property_name!r} -> source_feature={source_feature.name!r} "
                f"uid={source_feature.uid!r}"
            )
            matches.append(
                (property_name, source_feature, synced_property_name.strip())
            )
        if len(matches) == 1:
            logger.important(
                "notion sync backward-debug: infer success "
                f"feature={matches[0][0]!r} source_uid={matches[0][1].uid!r}"
            )
            return matches[0][0], matches[0][1]
        if len(matches) > 1:
            if sys.stdin is None or not sys.stdin.isatty():
                logger.warning(
                    "notion sync backward relation: multiple candidates found in "
                    "a non-interactive session; skipping auto-configuration"
                )
                logger.important(
                    "notion sync backward-debug: infer ambiguous candidates "
                    f"{[(name, feature.uid) for name, feature, _ in matches]}"
                )
                return None, None
            RICH_CONSOLE.print(
                "[bold yellow]notion backward relation[/] multiple candidates found; "
                "choose the feature to configure as backward-derived:",
                markup=True,
                highlight=False,
            )
            for i, (property_name, source_feature, synced_property_name) in enumerate(
                matches, start=1
            ):
                RICH_CONSOLE.print(
                    f"  {i}. {property_name} <-- {synced_property_name} "
                    f"(source uid: {source_feature.uid})",
                    markup=True,
                    highlight=False,
                )
            while True:
                selected = input(
                    "Choose backward relation candidate [1-"
                    f"{len(matches)}] (Enter to skip): "
                ).strip()
                if selected == "":
                    logger.important(
                        "notion sync backward-debug: user skipped ambiguous candidate selection"
                    )
                    return None, None
                if selected.isdigit():
                    idx = int(selected)
                    if 1 <= idx <= len(matches):
                        property_name, source_feature, _ = matches[idx - 1]
                        logger.important(
                            "notion sync backward-debug: user selected candidate "
                            f"feature={property_name!r} source_uid={source_feature.uid!r}"
                        )
                        return property_name, source_feature
                RICH_CONSOLE.print(
                    f"[yellow]Invalid choice. Use 1-{len(matches)} or press Enter to skip.[/]",
                    markup=True,
                    highlight=False,
                )
        logger.important(
            "notion sync backward-debug: infer no-unique-match "
            f"matches={[(name, feature.uid) for name, feature, _ in matches]}"
        )
        return None, None

    @staticmethod
    def _relation_feature_matches_target_type(
        local_feature: Any, target_type: Any
    ) -> bool:
        """Whether a local relation feature points to a specific target record type."""
        from lamindb.models.feature import parse_dtype

        dtype_str = getattr(local_feature, "_dtype_str", None)
        if not isinstance(dtype_str, str) or not dtype_str:
            return True
        parsed = parse_dtype(dtype_str)
        if len(parsed) != 1:
            return False
        parsed_dtype = parsed[0]
        if parsed_dtype.get("registry_str") != "Record":
            return False
        registry = parsed_dtype.get("registry")
        registry_uid = getattr(registry, "uid", None)
        target_uid = getattr(target_type, "uid", None)
        if isinstance(registry_uid, str) and isinstance(target_uid, str):
            return registry_uid == target_uid
        registry_name = getattr(registry, "name", None)
        target_name = getattr(target_type, "name", None)
        if isinstance(registry_name, str) and isinstance(target_name, str):
            return registry_name.lower() == target_name.lower()
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
        if schema is not None:
            for name, dtype_label, dtype in feature_plan:
                existing_feature = schema_members_by_name.get(name)
                if existing_feature is None:
                    continue
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

        if type_update_specs:
            for name, dtype_label, dtype, _ in type_update_specs:
                detail = self._feature_plan_detail(
                    db_name,
                    name,
                    dtype_label,
                    dtype,
                    record_field_mappings,
                    index_feature_name=index_feature_name,
                )
                if apply:
                    self._append_unique(report.updated_features, detail)
                else:
                    self._append_unique(report.update_features, detail)

        if apply and feature_type is not None and missing_specs:
            for name, _, dtype in missing_specs:
                ln.Feature(name=name, dtype=dtype, type=feature_type).save()
            logger.important(
                f"notion sync metadata: created {len(missing_specs)} features for {db_name!r}"
            )

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
            logger.important(
                f"notion sync metadata: updated {len(type_update_specs)} features to type {db_name!r}"
            )

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
        backward_feature_name: str | None = None
        backward_source_feature: Any | None = None
        if schema_spec:
            (
                backward_feature_name,
                backward_source_feature,
            ) = self._infer_notion_backward_relation_feature(
                schema_spec, features_by_name
            )
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
                        backward_feature_name is not None
                        and feature.name == backward_feature_name
                        and backward_source_feature is not None
                    ):
                        schema_features.append(
                            feature.with_config(
                                field=mapped_field, backward=backward_source_feature
                            )
                        )
                    elif mapped_field is None:
                        schema_features.append(feature)
                    else:
                        schema_features.append(feature.with_config(field=mapped_field))
                schema = ln.Schema(
                    schema_features,
                    name=db_name,
                    index=index_feature,
                ).save()
                self._append_unique(report.created_schemas, db_name)
                logger.important(f"notion sync metadata: created schema {db_name!r}")
            else:
                self._append_unique(report.create_schemas, db_name)
        else:
            logger.important(f"notion sync metadata: schema {db_name!r} already exists")
            if missing_specs:
                if apply:
                    self._append_unique(report.updated_schemas, db_name)
                else:
                    self._append_unique(report.update_schemas, db_name)
            existing_backward_uids = dict(schema._backward_feature_uids)
            target_feature = (
                features_by_name.get(backward_feature_name)
                if backward_feature_name is not None
                else None
            )
            backward_mapping_needed = (
                backward_source_feature is not None
                and target_feature is not None
                and existing_backward_uids.get(target_feature.uid)
                != backward_source_feature.uid
            )
            logger.important(
                "notion sync backward-debug: schema backward-eval "
                f"db={db_name!r} inferred_feature={backward_feature_name!r} "
                f"inferred_source_uid={getattr(backward_source_feature, 'uid', None)!r} "
                f"existing_backward_uids={existing_backward_uids!r} "
                f"target_feature_uid={getattr(target_feature, 'uid', None)!r} "
                f"needed={backward_mapping_needed}"
            )
            if backward_mapping_needed:
                if apply:
                    self._append_unique(report.updated_schemas, db_name)
                else:
                    self._append_unique(report.update_schemas, db_name)
            if apply and record_field_mappings:
                schema_record_fields = dict(schema._record_fields)
                changed = False
                for feature in features:
                    mapped_field = record_field_mappings.get(feature.name)
                    if mapped_field is None:
                        continue
                    if schema_record_fields.get(feature.uid) != mapped_field:
                        schema_record_fields[feature.uid] = mapped_field
                        changed = True
                if changed:
                    schema._aux = schema._aux or {}
                    schema._aux.setdefault("af", {})["2"] = schema_record_fields
                    schema.save(update_fields=["_aux"])
                    logger.important(
                        f"notion sync metadata: updated record-field mappings for schema {db_name!r}"
                    )
            if apply and backward_mapping_needed and target_feature is not None:
                backward_mappings = dict(schema._backward_feature_uids)
                backward_mappings[target_feature.uid] = backward_source_feature.uid
                schema._backward_feature_uids = backward_mappings
                schema.save(update_fields=["_aux"])
                logger.important(
                    "notion sync metadata: updated backward relation mapping "
                    f"for schema {db_name!r}"
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
                "notion sync backward-debug: loaded schema spec "
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
            "notion sync backward-debug: loaded schema spec "
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
        limit: int | None = None,
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
        db_id_set, parent_pages = self._collect_database_ids(parent_ids, limit=limit)
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
            if limit == 0:
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
            rec_type = self._resolve_record_type(
                db_id,
                apply=apply,
                report=report,
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
                rows = self.reader.rows(db_id, limit=limit, include_page_emoji=apply)
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
                    planned_transfers_by_db[db_id] = transfer_map
                    for detail in transfer_details:
                        if detail not in report.create_artifacts:
                            report.create_artifacts.append(detail)
                    continue

                before = _existing_by_ref(rec_type)
                before_edit = self._existing_edit_map(before)

                writes: list[dict[str, Any]] = []
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
                    else:
                        report.updated += 1
                        writes.append(row)
                to_write[db_id] = writes
                logger.important(
                    f"notion sync phase A: db={_compact_uuid(db_id)}, create={sum(1 for row in writes if row['notion_id'] not in before)}, "
                    f"update={sum(1 for row in writes if row['notion_id'] in before)}, unchanged={len(rows) - len(writes)}"
                )
                _, _, file_props = _kinds(db_specs[db_id])
                transfer_map, transfer_details = _planned_missing_file_transfers(
                    writes, file_props
                )
                planned_transfers_by_db[db_id] = transfer_map
                if not apply:
                    for detail in transfer_details:
                        if detail not in report.create_artifacts:
                            report.create_artifacts.append(detail)

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


@ln.flow("Ofbk5ruuTiN2")
def sync_from_notion(
    *,
    parents: str | list[str],
    token: str | None = None,
    apply: bool = False,
    limit: int | None = None,
) -> SyncReport:
    """Sync Notion pages via the class-based sync API."""
    syncer = _NotionSyncer(token=token)
    if isinstance(parents, str):
        parent_list = [parents]
    elif isinstance(parents, list):
        parent_list = list(parents)
    else:
        raise TypeError("parents must be a str or list[str].")
    report = syncer.import_pages(parents=parent_list, apply=apply, limit=limit)
    RICH_CONSOLE.print(report.to_pretty_text(), markup=True, highlight=False)
    return report


__all__ = [
    "sync_from_notion",
    "SyncReport",
]
