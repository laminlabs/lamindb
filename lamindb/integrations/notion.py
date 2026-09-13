"""Sync Notion pages to LaminDB records.

.. autofunction:: sync_from_notion
.. autoclass:: SyncReport

"""

from __future__ import annotations

import json
import os
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

import requests
from lamin_utils import logger

import lamindb as ln

API_VERSION = "2026-03-11"
BASE = "https://api.notion.com/v1"


@dataclass
class SyncReport:
    discovered: int = 0
    created: int = 0
    updated: int = 0
    unchanged: int = 0
    pending_relations: int = 0
    failed: int = 0
    errors: list[str] = field(default_factory=list)
    databases: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "discovered": self.discovered,
            "created": self.created,
            "updated": self.updated,
            "unchanged": self.unchanged,
            "pending_relations": self.pending_relations,
            "failed": self.failed,
            "errors": self.errors,
            "databases": self.databases,
        }


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


class _NotionReader:
    """Read-only Notion reader. Databases contain data sources; rows live on the data source."""

    def __init__(self, token: str) -> None:
        if not token:
            raise ValueError("A Notion access token is required.")
        self.s = requests.Session()
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
                raise LookupError(f"No data sources on database {database_id!r}.")
            if len(sources) > 1:
                logger.warning(
                    f"database {database_id!r} has {len(sources)} data sources; "
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
            row: dict[str, Any] = {"notion_id": None, "last_edited_time": None}
            for name, prop in page.get("properties", {}).items():
                if name in drop or prop.get("type") in drop:
                    continue
                row[name] = _flatten(prop)
            # page-level fields win over any same-named user property
            row["notion_id"] = page.get("id")
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
            ``last_edited_time``. Relation/people values are lists of Notion
            page UUIDs — the stable join key.
        """
        drop = drop or set()
        raw = self._call("GET", f"/pages/{page_id}")
        row: dict[str, Any] = {"notion_id": None, "last_edited_time": None}
        for name, prop in raw.get("properties", {}).items():
            if name in drop or prop.get("type") in drop:
                continue
            row[name] = _flatten(prop)
        row["notion_id"] = raw.get("id")
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
        if prop in ("notion_id", "last_edited_time") or val in (None, [], ""):
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

    if "notion_last_edited" in feat:
        values[feat["notion_last_edited"]] = row.get("last_edited_time")
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
        rec = by_id.get(nid)
        if rec is None:
            by_id[nid] = ln.Record(
                name=name or None,
                type=rec_type,
                reference=nid,
                reference_type="notion",
            ).save()
        elif name and rec.name != name:
            rec.name = name
            rec.save()
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

    def _collect_database_ids(self, parents: list[str]) -> set[str]:
        database_ids: set[str] = set()
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
                    f"Parent {parent!r} is neither a readable database nor page."
                )
            self._collect_databases_from_block(parent, seen_blocks, database_ids)
        return database_ids

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

    def _resolve_record_type(self, database_id: str):
        payload = self.reader._call("GET", f"/databases/{database_id}")
        db_name = self._database_title(payload, fallback=database_id)
        qs = ln.Record.filter(name=db_name, is_type=True)
        count = qs.count()
        if count == 0:
            raise ValueError(
                f"No Lamin record type named {db_name!r} for Notion database {database_id!r}."
            )
        if count > 1:
            raise ValueError(
                f"Ambiguous Lamin record type name {db_name!r}: found {count} matches."
            )
        return qs.one()

    def _validate_schema(self, database_id: str, rec_type) -> None:
        notion_props = set(self.reader.columns(database_id))
        schema_features = self._schema_feature_names(rec_type)
        missing_features = sorted(notion_props - schema_features)
        extra_features = sorted(schema_features - notion_props - {"notion_last_edited"})
        if "notion_last_edited" not in schema_features:
            raise ValueError(
                f"Record type {rec_type.name!r} is missing required feature "
                "'notion_last_edited'."
            )
        if missing_features or extra_features:
            problems: list[str] = []
            if missing_features:
                problems.append(f"missing in Lamin schema: {missing_features}")
            if extra_features:
                problems.append(f"extra in Lamin schema: {extra_features}")
            msg = "; ".join(problems)
            raise ValueError(
                f"Schema mismatch for database {database_id!r} <-> record type "
                f"{rec_type.name!r}: {msg}"
            )

    @staticmethod
    def _existing_edit_map(by_id: dict) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for notion_id, record in by_id.items():
            out[notion_id] = record.features.get_values().get("notion_last_edited")
        return out

    def import_pages(
        self,
        parents: str | list[str],
        *,
        dry_run: bool = False,
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

        report = SyncReport()
        db_ids = sorted(self._collect_database_ids(parent_ids))
        if not db_ids:
            raise ValueError(
                "No child databases discovered under parents. In phase 1, sync operates "
                "on page trees that include at least one Notion database."
            )
        report.databases = db_ids

        # Step 1: resolve and validate schema parity before any write.
        rec_types: dict[str, Any] = {}
        for db_id in db_ids:
            rec_type = self._resolve_record_type(db_id)
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

                before = _existing_by_ref(rec_type)
                before_edit = self._existing_edit_map(before)

                writes: list[dict[str, Any]] = []
                for row in rows:
                    notion_id = row["notion_id"]
                    edited = row.get("last_edited_time")
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

                if not dry_run:
                    after_maps[db_id] = _upsert_all(rec_type, rows)
                else:
                    after_maps[db_id] = before

            if dry_run:
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
    parents: list[str] | tuple[str, ...] | str,
    token: str | None = None,
    dry_run: bool = False,
    limit: int | None = None,
) -> SyncReport:
    """Sync Notion pages via the class-based sync API."""
    syncer = _NotionSyncer(token=token)
    if isinstance(parents, str):
        parent_list = [parents]
    else:
        parent_list = list(parents)
    report = syncer.import_pages(parents=parent_list, dry_run=dry_run, limit=limit)
    logger.important(f"{json.dumps(report.as_dict(), sort_keys=True)}")
    return report


__all__ = [
    "sync_from_notion",
    "SyncReport",
]
