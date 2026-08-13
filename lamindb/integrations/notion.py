"""Read a Notion database into Lamin records.

Two halves:

* :class:`Reader` — a read-only Notion client (databases -> data sources ->
  pages), flattening pages to plain dicts. Relation values are always lists of
  Notion page UUIDs (the stable join key), never titles.
* the sync layer — :func:`sync_all` pulls several databases in the efficient
  order (upsert everything, then write everything, so every cross-database link
  resolves in a single pass). :func:`import_db` / :func:`link` are the per-
  database path for incremental updates afterwards.
"""

from __future__ import annotations

import json
import time
from contextlib import contextmanager
from typing import Any

import requests
from lamin_utils import logger

API_VERSION = "2026-03-11"
BASE = "https://api.notion.com/v1"


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


class Reader:
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
    import lamindb as ln

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
    import lamindb as ln

    return {
        r.reference: r for r in ln.Record.filter(type=rec_type, reference_type="notion")
    }


def _resolved_map(uuids) -> dict:
    """{notion_uuid: ln.Record} for the UUIDs that resolve — ONE query, not N."""
    import lamindb as ln

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
    import lamindb as ln

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


def upsert(rec_type, notion_id: str, name: str | None = None):
    """Create-or-get a single record keyed on its Notion UUID (stable identity).

    A changed title is a plain attribute update on the already-saved row (an
    UPDATE, not a re-construction), so it cannot trigger name-based merge. For
    bulk loads prefer :func:`import_db` / :func:`sync_all`, which resolve the
    existing-record map once instead of querying per row.
    """
    import lamindb as ln

    rec = ln.Record.filter(
        type=rec_type, reference=notion_id, reference_type="notion"
    ).one_or_none()
    if rec is None:
        return ln.Record(
            name=name or None,
            type=rec_type,
            reference=notion_id,
            reference_type="notion",
        ).save()
    if name and rec.name != name:
        rec.name = name
        rec.save()
    return rec


def materialize(record, row: dict, spec: dict, resolved: dict, prop_map=None) -> int:
    """Rebuild ONE record's full feature state from its Notion row; write once.

    Standalone helper — resolves the schema and label kinds itself. Bulk paths
    use the internal fast loop instead (schema resolved once for all rows).
    Returns the count of relation targets that did not resolve (still pending).
    """
    rel, lab = _kinds(spec)
    schema = (
        record.type.schema
        if (record.type is not None and record.type.schema is not None)
        else record.schema
    )
    feat = _feat_map(schema)
    values, pending = _row_values(
        row, rel, lab, feat, resolved, prop_map, create_labels=True
    )
    record.features.set_values(values)
    return pending


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
    import lamindb as ln

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


def import_db(reader, database_id: str, rec_type, drop=None, prop_map=None) -> dict:
    """Import one Notion database: upsert every row, then write it once.

    Same-database relations (e.g. people->people) resolve immediately. Relations
    into a database you haven't imported yet come back as ``pending`` — import
    that database and call :func:`link`, or use :func:`sync_all` to do the whole
    set in one efficient pass. Returns ``{"records": int, "pending": int}``.
    """
    with _bulk_creation():
        spec = reader.schema(database_id)
        rows = reader.rows(database_id, drop=drop)  # single paginated pull
        by_id = _upsert_all(rec_type, rows)
        return _write(reader, rows, rec_type, spec, prop_map, by_id=by_id)


def link(reader, database_id: str, rec_type, drop=None, prop_map=None) -> dict:
    """Re-resolve one database's relations against everything imported so far.

    Idempotent and convergent: run it after importing more databases and
    previously-``pending`` links attach. Re-pulls rows because a full-replace
    write needs the scalar values too. Returns ``{"records", "pending"}``.
    """
    with _bulk_creation():
        spec = reader.schema(database_id)
        rows = reader.rows(database_id, drop=drop)
        return _write(reader, rows, rec_type, spec, prop_map)


def sync_all(reader, mapping: dict, drop=None, prop_maps=None) -> dict:
    """Import several databases in the efficient order — no re-link, no double write.

    Phase A upserts every row of every database, so all references exist. Phase B
    then writes each database exactly once, and because every target already
    exists, cross-database links (people->org, meeting->people) resolve in that
    single pass. This replaces the import-then-relink dance and halves the writes.

    Args:
        mapping: ``{rec_type: database_id}`` — the Record type for each database.
        drop: property names/types to omit on every page.
        prop_maps: ``{rec_type: {notion_prop: feature_name}}`` per database.

    Returns:
        ``{rec_type: {"records": int, "pending": int}}``. Any leftover ``pending``
        means a relation target lives outside the databases you passed.
    """
    prop_maps = prop_maps or {}
    with _bulk_creation():
        specs: dict = {}
        all_rows: dict = {}
        maps: dict = {}
        # Phase A — pull + upsert every database so all references exist
        for rec_type, db_id in mapping.items():
            specs[rec_type] = reader.schema(db_id)
            all_rows[rec_type] = reader.rows(db_id, drop=drop)
            maps[rec_type] = _upsert_all(rec_type, all_rows[rec_type])
        # Phase B — write each database once; every link now resolves
        out: dict = {}
        for rec_type, _db_id in mapping.items():
            out[rec_type] = _write(
                reader,
                all_rows[rec_type],
                rec_type,
                specs[rec_type],
                prop_maps.get(rec_type),
                by_id=maps[rec_type],
            )
        return out


__all__ = ["Reader", "upsert", "materialize", "import_db", "link", "sync_all"]
