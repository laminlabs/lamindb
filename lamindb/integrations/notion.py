"""Read a Notion database into a JSON table."""

from __future__ import annotations

import json
from typing import Any

import requests

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

    def _call(self, method: str, path: str, body: dict | None = None) -> dict:
        r = self.s.request(method, f"{BASE}{path}", json=body, timeout=30)
        if r.status_code == 401:
            raise PermissionError("Invalid or expired Notion token.")
        if r.status_code == 404:
            raise LookupError(
                f"404 on {path} — not found, or not shared with this connection "
                "(Notion: ••• -> Connections -> Connect to)."
            )
        r.raise_for_status()
        return r.json()

    def data_sources(self, database_id: str) -> list[dict]:
        """List the data sources under a database: [{id, name}, ...]."""
        return self._call("GET", f"/databases/{database_id}").get("data_sources", [])

    def _resolve(self, database_id: str) -> str:
        if database_id not in self._ds:
            sources = self.data_sources(database_id)
            if not sources:
                raise LookupError(f"No data sources on database {database_id!r}.")
            if len(sources) > 1:
                print(
                    f"note: {len(sources)} data sources; using {sources[0]['name']!r}"
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

    def _query(self, ds: str) -> list[dict]:
        """Every page in a data source, raw. Paginates until exhausted."""
        body: dict[str, Any] = {"page_size": 100}
        pages: list[dict] = []
        while True:
            payload = self._call("POST", f"/data_sources/{ds}/query", body)
            pages.extend(payload.get("results", []))
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
                print(f"  ! data source {ds_id} not shared — leaving IDs in place")
                self._titles[ds_id] = {}
        return self._titles[ds_id]

    def rows(
        self,
        database_id: str,
        drop: set[str] | None = None,
    ) -> list[dict]:
        """Every page flattened to a dict.

        Relation values are always lists of Notion page UUIDs — never titles.
        UUIDs are the stable join key; titles are not unique and change on
        rename. Use :meth:`relation_titles` for display names.

        Args:
            database_id: the Notion database ID.
            drop: property names or property types to omit, e.g.
                {"messages", "created_by"}.
        """
        drop = drop or set()
        ds = self._resolve(database_id)

        rows: list[dict] = []
        for page in self._query(ds):
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
