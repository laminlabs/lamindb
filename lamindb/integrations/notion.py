"""Sync Notion pages to LaminDB records.

.. autofunction:: sync_from_notion
.. autoclass:: SyncReport

"""

from __future__ import annotations

import json
import os
import re
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import Any

import httpx
from lamin_utils import logger
from rich.console import Console

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
    created_feature_types: list[str] = field(default_factory=list)
    create_feature_types: list[str] = field(default_factory=list)
    created_schemas: list[str] = field(default_factory=list)
    create_schemas: list[str] = field(default_factory=list)
    created_features: list[str] = field(default_factory=list)
    create_features: list[str] = field(default_factory=list)

    def to_pretty_text(self) -> str:
        """Render a concise human-readable sync report."""

        def metric(key: str, value: str | int, color: str = "white") -> str:
            return f"[bold white]{key}[/]: [{color}]{value}[/]"

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
        if self.create_schemas:
            lines.append(
                metric("create_schemas", ", ".join(self.create_schemas), action_color)
            )
        if self.created_schemas:
            lines.append(
                metric("created_schemas", ", ".join(self.created_schemas), "green")
            )
        if self.create_features:
            lines.append("[bold]create_features[/]:")
            lines.extend(
                f"  [{action_color}]{feature}[/]" for feature in self.create_features
            )
        if self.created_features:
            lines.append("[bold]created_features[/]:")
            lines.extend(f"  [green]{feature}[/]" for feature in self.created_features)
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
        """{property_name: {"type": str, "target": str | None, "dual": dict | None}}.

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
            if t == "relation":
                rel = p.get("relation", {})
                target = rel.get("data_source_id") or rel.get("database_id")
                dual = rel.get("dual_property")
            out[name] = {"type": t, "target": target, "dual": dual}
        self._schema[database_id] = out
        return out

    def columns(self, database_id: str) -> dict[str, str]:
        """Return {property_name: notion_type}."""
        return {k: v["type"] for k, v in self.schema(database_id).items()}

    def _query(self, ds: str, limit: int | None = None) -> list[dict]:
        """Pages in a data source, raw. Paginates until exhausted or `limit` reached."""
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
            rows.append(row)
        return rows

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


def _kinds(spec: dict) -> tuple[set, set]:
    """Split a Notion schema into (relation-props, label-props)."""
    rel = {p for p, s in spec.items() if s["type"] in ("relation", "people")}
    lab = {
        p for p, s in spec.items() if s["type"] in ("select", "status", "multi_select")
    }
    return rel, lab


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


def _row_values(row, rel, lab, feat, resolved, prop_map, create_labels):
    """Build the full {Feature: value} dict for one row. Returns (values, pending)."""
    prop_map = prop_map or {}
    values: dict[Any, Any] = {}
    pending = 0
    for prop, val in row.items():
        if prop in ("notion_id", "created_time", "last_edited_time") or val in (
            None,
            [],
            "",
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
        else:
            values[f] = val
    return values, pending


def _write(reader, rows, rec_type, spec, prop_map=None, by_id=None) -> dict:
    """Materialize every row of one database. Schema, kinds and labels resolved once."""
    rel, lab = _kinds(spec)
    feat = _feat_map(rec_type.schema)

    _batch_labels(rows, lab)  # every ULabel created in one call

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
            row, rel, lab, feat, resolved, prop_map, create_labels=False
        )
        rec.features.set_values(values)
        pending += p
        records += 1
    return {"records": records, "pending": pending}


def _upsert_all(rec_type, rows) -> dict:
    """Upsert every row of one database against a single existing-record map.

    Returns the {notion_uuid: record} map (existing + newly created), ready to
    hand to :func:`_write` so it never re-queries.
    """
    by_id = _existing_by_ref(rec_type)  # ONE query, not one per row
    for row in rows:
        nid, name = row["notion_id"], row.get("name")
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
        self, block_id: str, seen: set[str], out: set[str]
    ) -> None:
        for block in self._iter_block_children(block_id):
            bid = block.get("id")
            if bid and bid in seen:
                continue
            if bid:
                seen.add(bid)
            if block.get("type") == "child_database" and bid:
                out.add(bid)
            if block.get("has_children") and bid:
                self._collect_databases_from_block(bid, seen, out)

    def _collect_database_ids(
        self, parents: list[str]
    ) -> tuple[set[str], dict[str, str]]:
        database_ids: set[str] = set()
        parent_pages: dict[str, str] = {}
        seen_blocks: set[str] = set()
        for parent in parents:
            db_payload = self._safe_call(f"/databases/{parent}")
            if db_payload is not None:
                database_ids.add(parent)
                # recurse through rows as pages to discover nested child databases
                for row in self.reader.rows(parent):
                    notion_id = row.get("notion_id")
                    if notion_id:
                        self._collect_databases_from_block(
                            notion_id, seen_blocks, database_ids
                        )
                continue
            page_payload = self._safe_call(f"/pages/{parent}")
            if page_payload is None:
                raise LookupError(
                    f"Parent {_compact_uuid(parent)!r} is neither a readable database nor page."
                )
            parent_title = _page_title(page_payload).strip() or _compact_uuid(parent)
            parent_pages[parent] = parent_title
            self._collect_databases_from_block(parent, seen_blocks, database_ids)
        return database_ids, parent_pages

    def _resolve_or_create_type_by_name(
        self, name: str, *, apply: bool, report: SyncReport
    ):
        qs = ln.Record.filter(name=name, is_type=True)
        count = qs.count()
        if count == 0:
            if apply:
                ln.Record(name=name, is_type=True).save()
                if name not in report.created_record_types:
                    report.created_record_types.append(name)
            elif name not in report.create_record_types:
                report.create_record_types.append(name)
            return
        if count > 1:
            raise ValueError(
                f"Ambiguous Lamin record type name {name!r}: found {count} matches."
            )

    @staticmethod
    def _database_title(payload: dict, fallback: str) -> str:
        title = payload.get("title") or []
        text = "".join(part.get("plain_text", "") for part in title).strip()
        return text or fallback

    @staticmethod
    def _schema_feature_names(rec_type) -> set[str]:
        schema = rec_type.schema
        if schema is None:
            raise ValueError(
                f"Record type {rec_type.name!r} has no schema. Add a schema before syncing."
            )
        return {feature.name for feature in schema.members}

    @staticmethod
    def _feature_dtype_from_notion_type(notion_type: str):
        if notion_type == "number":
            return "num"
        if notion_type == "checkbox":
            return bool
        if notion_type in {"multi_select", "people", "relation", "files"}:
            return list[str]
        return str

    @staticmethod
    def _feature_dtype_label_from_notion_type(notion_type: str) -> str:
        if notion_type == "number":
            return "num"
        if notion_type == "checkbox":
            return "bool"
        if notion_type in {"multi_select", "people", "relation", "files"}:
            return "list[str]"
        return "str"

    @staticmethod
    def _index_feature_name_from_columns(columns: dict[str, str]) -> str | None:
        for name, notion_type in columns.items():
            if notion_type == "title":
                return name
        return None

    def _database_feature_plan(
        self, database_id: str, columns: dict[str, str] | None = None
    ) -> list[tuple[str, str, Any]]:
        if columns is None:
            columns = self.reader.columns(database_id)
        ordered_feature_names = list(columns)
        plan: list[tuple[str, str, Any]] = []
        for name in ordered_feature_names:
            notion_type = columns[name]
            plan.append(
                (
                    name,
                    self._feature_dtype_label_from_notion_type(notion_type),
                    self._feature_dtype_from_notion_type(notion_type),
                )
            )
        return plan

    @staticmethod
    def _append_unique(values: list[str], value: str) -> None:
        if value not in values:
            values.append(value)

    def _plan_or_create_db_metadata(
        self,
        db_name: str,
        feature_plan: list[tuple[str, str, Any]],
        index_feature_name: str | None = None,
        *,
        apply: bool,
        report: SyncReport,
    ) -> tuple[Any, list[Any], Any]:
        feature_type_qs = ln.Feature.filter(name=db_name, is_type=True)
        feature_type_count = feature_type_qs.count()
        if feature_type_count > 1:
            raise ValueError(
                f"Ambiguous LaminDB feature type name {db_name!r}: found {feature_type_count} matches."
            )
        feature_type = feature_type_qs.one_or_none()
        if feature_type is None:
            if apply:
                feature_type = ln.Feature(name=db_name, is_type=True).save()
                self._append_unique(report.created_feature_types, db_name)
            else:
                self._append_unique(report.create_feature_types, db_name)

        missing_specs: list[tuple[str, str, Any]] = []
        if feature_type is None:
            missing_specs = feature_plan
        else:
            existing_names = set(
                ln.Feature.filter(
                    name__in=[name for name, _, _ in feature_plan], type=feature_type
                ).values_list("name", flat=True)
            )
            missing_specs = [
                spec for spec in feature_plan if spec[0] not in existing_names
            ]

        if missing_specs:
            for name, dtype_label, _ in missing_specs:
                detail = f"{db_name} / {name}: {dtype_label}"
                if apply:
                    self._append_unique(report.created_features, detail)
                else:
                    self._append_unique(report.create_features, detail)

        if apply and feature_type is not None and missing_specs:
            for name, _, dtype in missing_specs:
                ln.Feature(name=name, dtype=dtype, type=feature_type).save()

        if feature_type is not None:
            features = list(
                ln.Feature.filter(
                    name__in=[name for name, _, _ in feature_plan], type=feature_type
                )
            )
        else:
            features = []

        schema_qs = ln.Schema.filter(name=db_name)
        schema_count = schema_qs.count()
        if schema_count > 1:
            raise ValueError(
                f"Ambiguous LaminDB schema name {db_name!r}: found {schema_count} matches."
            )
        schema = schema_qs.one_or_none()
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
                schema = ln.Schema(
                    features,
                    name=db_name,
                    index=index_feature,
                ).save()
                self._append_unique(report.created_schemas, db_name)
            else:
                self._append_unique(report.create_schemas, db_name)
        return feature_type, features, schema

    def _create_record_type(self, database_id: str, db_name: str, report: SyncReport):
        columns = self.reader.columns(database_id)
        feature_plan = self._database_feature_plan(database_id, columns=columns)
        index_feature_name = self._index_feature_name_from_columns(columns)
        _, _, schema = self._plan_or_create_db_metadata(
            db_name,
            feature_plan,
            index_feature_name=index_feature_name,
            apply=True,
            report=report,
        )
        assert schema is not None  # schema is always created/resolved in apply mode
        return ln.Record(name=db_name, is_type=True, schema=schema).save()

    def _resolve_record_type(
        self, database_id: str, *, apply: bool, report: SyncReport
    ):
        payload = self.reader._call("GET", f"/databases/{database_id}")
        db_name = self._database_title(payload, fallback=database_id)
        qs = ln.Record.filter(name=db_name, is_type=True)
        count = qs.count()
        if count == 0:
            columns = self.reader.columns(database_id)
            feature_plan = self._database_feature_plan(database_id, columns=columns)
            index_feature_name = self._index_feature_name_from_columns(columns)
            if not apply:
                report.create_record_types.append(db_name)
                self._plan_or_create_db_metadata(
                    db_name,
                    feature_plan,
                    index_feature_name=index_feature_name,
                    apply=False,
                    report=report,
                )
                return None
            rec_type = self._create_record_type(database_id, db_name, report=report)
            report.created_record_types.append(db_name)
            return rec_type
        if count > 1:
            raise ValueError(
                f"Ambiguous Lamin record type name {db_name!r}: found {count} matches."
            )
        return qs.one()

    def _validate_schema(self, database_id: str, rec_type) -> None:
        notion_props = set(self.reader.columns(database_id))
        schema_features = self._schema_feature_names(rec_type)
        missing_features = sorted(notion_props - schema_features)
        extra_features = sorted(schema_features - notion_props)
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
        db_id_set, parent_pages = self._collect_database_ids(parent_ids)
        db_ids = sorted(db_id_set)
        if not db_ids:
            raise ValueError(
                "No child databases discovered under parents. In phase 1, sync operates "
                "on page trees that include at least one Notion database."
            )
        report.databases = [_compact_uuid(db_id) for db_id in db_ids]

        # Parent pages can also map to LaminDB record types.
        for parent_type_name in sorted(set(parent_pages.values())):
            self._resolve_or_create_type_by_name(
                parent_type_name, apply=apply, report=report
            )

        # Step 1: resolve and validate schema parity before any write.
        rec_types: dict[str, Any] = {}
        for db_id in db_ids:
            rec_type = self._resolve_record_type(db_id, apply=apply, report=report)
            if rec_type is not None:
                self._validate_schema(db_id, rec_type)
            rec_types[db_id] = rec_type

        after_maps: dict[str, dict[str, Any]] = {}
        to_write: dict[str, list[dict[str, Any]]] = {}

        with _bulk_creation():
            # Phase A: discover + upsert identity rows.
            for db_id in db_ids:
                rec_type = rec_types[db_id]
                rows = self.reader.rows(db_id, limit=limit)
                report.discovered += len(rows)
                if rec_type is None:
                    # dry-run mode with missing type: all discovered rows are new.
                    report.created += len(rows)
                    to_write[db_id] = []
                    after_maps[db_id] = {}
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

                if apply:
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
                spec = self.reader.schema(db_id)
                write_rows = to_write[db_id]
                if not write_rows:
                    continue
                stats = _write(
                    self.reader,
                    write_rows,
                    rec_type,
                    spec,
                    prop_map=None,
                    by_id=after_maps[db_id],
                )
                report.pending_relations += stats["pending"]

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
