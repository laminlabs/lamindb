from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd
    import requests

# Pinned to 2022-06-28. Notion's 2025-09-03 version replaces the
# /v1/databases query endpoints with /v1/data_sources, since a database can
# now contain multiple data sources with distinct schemas. This client targets
# single-data-source databases; against a multi-source database the
# /v1/databases endpoints fail with a validation error. Migration path:
# GET /v1/databases/{id} -> data_sources[0].id -> POST /v1/data_sources/{id}/query.
NOTION_API_VERSION = "2022-06-28"
_NOTION_BASE_URL = "https://api.notion.com/v1"
_MAX_RETRIES = 3
_TIMEOUT = 30
_MAX_BACKOFF = 60
_RETRY_STATUSES = frozenset({429, 502, 503, 504})

_log = logging.getLogger(__name__)


def _flatten_property(prop: dict) -> Any:
    """Flatten a single Notion page property value to a plain Python object.

    Dispatches on ``prop["type"]`` and returns scalars, lists, or ``None``
    for derived / unrecognized types.

    For ``date`` properties only the ``start`` value is returned; ``end`` and
    ``time_zone`` are silently dropped.

    Args:
        prop: A raw Notion property value dict as returned by the Pages API.

    Returns:
        A plain Python scalar or list, or ``None`` for derived/unknown types.
    """
    ptype = prop.get("type")
    if ptype in ("title", "rich_text"):
        return "".join(seg.get("plain_text", "") for seg in (prop.get(ptype) or []))
    if ptype in ("email", "phone_number", "url", "created_time", "last_edited_time"):
        return prop.get(ptype)
    if ptype in ("number", "checkbox"):
        return prop.get(ptype)
    if ptype in ("select", "status"):
        option = prop.get(ptype)
        return option["name"] if option else None
    if ptype == "multi_select":
        return [opt["name"] for opt in (prop.get("multi_select") or [])]
    if ptype == "date":
        date_val = prop.get("date")
        return date_val["start"] if date_val else None
    if ptype == "people":
        return [person["id"] for person in (prop.get("people") or [])]
    if ptype == "relation":
        return [rel["id"] for rel in (prop.get("relation") or [])]
    if ptype == "files":
        # Notion's internal file URLs are signed and expire after ~1 hour.
        # Do not store them as durable pointers; re-fetch as needed.
        urls: list[str] = []
        for f in prop.get("files") or []:
            ftype = f.get("type")
            if ftype == "external":
                urls.append(f["external"]["url"])
            elif ftype == "file":
                urls.append(f["file"]["url"])
        return urls
    if ptype in ("rollup", "formula"):
        return None
    _log.debug("_flatten_property: unrecognized Notion property type %r", ptype)
    return None


class Notion:
    """Read-only client for the Notion API.

    Reads Notion databases and returns plain Python / pandas objects.
    No lamindb ORM imports; this module is intentionally registry-free and
    can be used standalone.

    Args:
        token: Notion integration secret token. Must be provided explicitly;
            no environment-variable fallback is applied.

    Raises:
        ValueError: If ``token`` is empty or falsy.

    Example:

        ::

            from lamindb.integrations import Notion

            df = Notion(token="<token>").to_dataframe("<database_id>")
    """

    def __init__(self, token: str) -> None:
        import requests as _requests

        if not token:
            raise ValueError("A Notion integration token is required.")
        self._session: requests.Session = _requests.Session()
        self._session.headers.update(
            {
                "Authorization": f"Bearer {token}",
                "Notion-Version": NOTION_API_VERSION,
                "Content-Type": "application/json",
            }
        )

    def _request(self, method: str, url: str, **kwargs: Any) -> requests.Response:
        # Retry is flat (no backoff growth) — sufficient for a read-only client.
        kwargs.setdefault("timeout", _TIMEOUT)
        for attempt in range(_MAX_RETRIES):
            resp = self._session.request(method, url, **kwargs)
            if resp.status_code in _RETRY_STATUSES and attempt < _MAX_RETRIES - 1:
                try:
                    retry_after = int(resp.headers.get("Retry-After", 5))
                except ValueError:
                    retry_after = 5
                retry_after = min(max(retry_after, 0), _MAX_BACKOFF)
                time.sleep(retry_after)
                continue
            break
        return resp

    def _raise_for_status(self, resp: requests.Response, action: str) -> None:
        if resp.ok:
            return
        if resp.status_code == 401:
            raise PermissionError(
                f"Notion API returned 401 while {action}: the token is invalid or "
                "expired. Check your Notion integration token."
            )
        if resp.status_code == 404:
            raise LookupError(
                f"Notion API returned 404 while {action}. "
                "The database may not be shared with the integration — "
                "open the Notion page, click Share, and invite the integration."
            )
        if resp.status_code == 429:
            raise RuntimeError(
                f"Notion API rate limit persisted while {action} after "
                f"{_MAX_RETRIES} attempts. Retry later or reduce concurrency."
            )
        resp.raise_for_status()

    def get_schema(self, database_id: str) -> dict[str, dict]:
        """Retrieve the schema of a Notion database.

        Calls ``GET /v1/databases/{id}``. Types come from the schema, not
        inferred by sampling pages.

        Args:
            database_id: The Notion database ID.

        Returns:
            A dict mapping each property name to a metadata dict with keys:

            - ``id`` (str): the Notion property ID (stable across renames)
            - ``type`` (str): the Notion property type
            - ``relation_database_id`` (str | None): linked database ID for
              ``relation`` properties, ``None`` otherwise
            - ``options`` (list[str] | None): option names for ``select``,
              ``multi_select``, and ``status`` properties, ``None`` otherwise

        Raises:
            LookupError: If the database is not found or not shared with the
                integration.
            PermissionError: If the token is invalid.
        """
        url = f"{_NOTION_BASE_URL}/databases/{database_id}"
        resp = self._request("GET", url)
        self._raise_for_status(resp, f"fetching schema for database {database_id!r}")
        schema: dict[str, dict] = {}
        for name, prop in resp.json().get("properties", {}).items():
            ptype = prop.get("type", "")
            relation_database_id: str | None = None
            options: list[str] | None = None
            if ptype == "relation":
                relation_database_id = prop.get("relation", {}).get("database_id")
            if ptype in ("select", "multi_select", "status"):
                raw_opts = prop.get(ptype, {}).get("options", [])
                options = [opt["name"] for opt in raw_opts]
            schema[name] = {
                "id": prop.get("id"),
                "type": ptype,
                "relation_database_id": relation_database_id,
                "options": options,
            }
        return schema

    def query_database(
        self,
        database_id: str,
        since: str | None = None,
    ) -> list[dict]:
        """Query all pages in a Notion database.

        Calls ``POST /v1/databases/{id}/query`` with ``page_size=100`` and
        follows ``has_more`` / ``next_cursor`` until exhausted.

        Args:
            database_id: The Notion database ID.
            since: ISO 8601 timestamp. When given, only pages with
                ``last_edited_time`` on or after this value are returned,
                enabling incremental syncs.

        Returns:
            Raw Notion page dicts, unmodified.

        Raises:
            LookupError: If the database is not found or not shared with the
                integration.
            PermissionError: If the token is invalid.
        """
        url = f"{_NOTION_BASE_URL}/databases/{database_id}/query"
        body: dict[str, Any] = {"page_size": 100}
        if since is not None:
            body["filter"] = {
                "timestamp": "last_edited_time",
                "last_edited_time": {"on_or_after": since},
            }
        pages: list[dict] = []
        while True:
            resp = self._request("POST", url, json=body)
            self._raise_for_status(resp, f"querying database {database_id!r}")
            payload = resp.json()
            pages.extend(payload.get("results", []))
            if not payload.get("has_more"):
                break
            body["start_cursor"] = payload["next_cursor"]
        return pages

    def to_dataframe(
        self,
        database_id: str,
        since: str | None = None,
    ) -> pd.DataFrame:
        """Fetch a Notion database and return it as a pandas DataFrame.

        Requires `pandas`, which ships with the ``full`` extra rather than the
        core install. :meth:`get_schema` and :meth:`query_database` have no
        such requirement.

        Each row is one Notion page. Property values are flattened to plain
        Python scalars or lists via ``_flatten_property``. Two columns are
        always the first two regardless of the database schema:

        - ``notion_id``: the page UUID (idempotency key for downstream sync)
        - ``last_edited_time``: page-level edit timestamp (incremental cursor)

        These two columns always reflect the page-level fields. If the
        database has a user property with the same name, the page-level value
        wins.

        Args:
            database_id: The Notion database ID.
            since: ISO 8601 timestamp forwarded to :meth:`query_database` to
                enable incremental syncs.

        Returns:
            A ``pd.DataFrame`` with one row per page. Returns an empty
            DataFrame with correct columns (in the same order as a populated
            result) when no pages match.

        Raises:
            LookupError: If the database is not found or not shared with the
                integration.
            PermissionError: If the token is invalid.
        """
        import pandas as pd

        pages = self.query_database(database_id, since=since)
        if not pages:
            schema_cols = list(self.get_schema(database_id))
            return pd.DataFrame(columns=["notion_id", "last_edited_time", *schema_cols])

        rows: list[dict[str, Any]] = []
        for page in pages:
            # Initialize with sentinel so notion_id and last_edited_time are
            # always the first two columns; reassignment below preserves position.
            row: dict[str, Any] = {"notion_id": None, "last_edited_time": None}
            for name, prop in page.get("properties", {}).items():
                row[name] = _flatten_property(prop)
            row["notion_id"] = page.get("id")
            row["last_edited_time"] = page.get("last_edited_time")
            rows.append(row)
        return pd.DataFrame(rows)


__all__ = ["Notion"]
