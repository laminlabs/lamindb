"""Tests for lamindb.integrations.notion.

Most tests mock HTTP at the httpx.Client level; one smoke test can run live
when NOTION_TOKEN is present.
"""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import httpx
import lamindb as ln
import pytest
from lamindb.integrations.notion import (
    API_VERSION,
    BASE,
    SyncReport,
    _artifact_key_from_url,
    _attach_page_markdown,
    _ensure_artifacts,
    _flatten,
    _NotionReader,
    _NotionSyncer,
    _planned_missing_file_transfers,
    _short_file_source,
    _upsert_all,
    _write,
    sync_from_notion,
)
from rich.console import Console

FIXTURES = Path(__file__).parent / "notion_test_data"


def _load_fixture(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text())


def _make_response(data: dict, status_code: int = 200) -> MagicMock:
    resp = MagicMock()
    resp.status_code = status_code
    resp.ok = status_code < 400
    resp.json.return_value = data
    return resp


DB = _load_fixture("notion_database.json")
DS = _load_fixture("notion_data_source.json")
PAGES = _load_fixture("notion_page.json")
DS_ID = "ds-1111-2222"

ORG_PAGES = {
    "results": [
        {
            "id": "related-page-id-111",
            "properties": {
                "Name": {"type": "title", "title": [{"plain_text": "Deepmind"}]}
            },
        }
    ],
    "has_more": False,
    "next_cursor": None,
}


@pytest.fixture()
def reader():
    """_NotionReader with httpx.Client replaced by a MagicMock."""
    with patch("httpx.Client") as MockSession:
        MockSession.return_value = MagicMock()
        client = _NotionReader(token="secret-test-token")  # noqa: S106
    return client


@pytest.fixture(scope="module")
def page_props():
    return PAGES["results"][0]["properties"]


# ---------------------------------------------------------------------------
# _flatten — one test per type in the dispatch table
# ---------------------------------------------------------------------------


def test_flatten_title(page_props):
    assert _flatten(page_props["Name"]) == "Row One"


def test_flatten_title_empty():
    assert _flatten({"type": "title", "title": []}) == ""


def test_flatten_rich_text(page_props):
    assert _flatten(page_props["Notes"]) == "some notes"


def test_flatten_rich_text_multipart():
    prop = {
        "type": "rich_text",
        "rich_text": [{"plain_text": "foo"}, {"plain_text": " bar"}],
    }
    assert _flatten(prop) == "foo bar"


def test_flatten_email(page_props):
    assert _flatten(page_props["Email"]) == "alice@example.com"


def test_flatten_email_null():
    assert _flatten({"type": "email", "email": None}) is None


def test_flatten_phone_number(page_props):
    assert _flatten(page_props["Phone"]) == "+1-555-0100"


def test_flatten_url(page_props):
    assert _flatten(page_props["Website"]) == "https://example.com"


def test_flatten_number(page_props):
    assert _flatten(page_props["Score"]) == 7.5


def test_flatten_checkbox_true(page_props):
    assert _flatten(page_props["Active"]) is True


def test_flatten_checkbox_false():
    assert _flatten({"type": "checkbox", "checkbox": False}) is False


def test_flatten_select(page_props):
    assert _flatten(page_props["Priority"]) == "High"


def test_flatten_select_null():
    assert _flatten({"type": "select", "select": None}) is None


def test_flatten_status(page_props):
    assert _flatten(page_props["State"]) == "Done"


def test_flatten_status_null():
    assert _flatten({"type": "status", "status": None}) is None


def test_flatten_multi_select(page_props):
    assert _flatten(page_props["Tags"]) == ["python", "data"]


def test_flatten_multi_select_empty():
    assert _flatten({"type": "multi_select", "multi_select": []}) == []


def test_flatten_date_returns_start_drops_end(page_props):
    assert _flatten(page_props["Due"]) == "2024-02-01"


def test_flatten_date_null():
    assert _flatten({"type": "date", "date": None}) is None


def test_flatten_people(page_props):
    assert _flatten(page_props["Owner"]) == ["user-id-abc"]


def test_flatten_relation_returns_page_ids(page_props):
    assert _flatten(page_props["Related"]) == ["related-page-id-111"]


def test_flatten_created_by(page_props):
    assert _flatten(page_props["Author"]) == "user-created-1"


def test_flatten_created_by_null():
    assert _flatten({"type": "created_by", "created_by": None}) is None


def test_flatten_last_edited_by(page_props):
    assert _flatten(page_props["LastEditor"]) == "user-edited-2"


def test_flatten_files(page_props):
    assert _flatten(page_props["Attachment"]) == [
        "https://example.com/doc.pdf",
        "https://s3.amazonaws.com/img.png",
    ]


def test_flatten_files_unknown_subtype_skipped():
    prop = {
        "type": "files",
        "files": [{"type": "future_cloud", "future_cloud": {"u": "x"}}],
    }
    assert _flatten(prop) == []


def test_planned_missing_file_transfers_lists_each_missing_file():
    rows = [
        {
            "notion_id": "row-1",
            "Attachment": ["https://example.com/a.pdf", "https://example.com/b.pdf"],
        }
    ]
    transfer_map, details = _planned_missing_file_transfers(rows, {"Attachment"})
    assert set(transfer_map) == {
        "https://example.com/a.pdf",
        "https://example.com/b.pdf",
    }
    assert details == [
        f"{_short_file_source('https://example.com/a.pdf')} <- row-1:Attachment",
        f"{_short_file_source('https://example.com/b.pdf')} <- row-1:Attachment",
    ]


def test_ensure_artifacts_downloads_then_saves_local_file():
    report = SyncReport(apply=True)
    created_artifact = type("ArtifactStub", (), {"uid": "abc"})()
    record = MagicMock()
    record.save.return_value = created_artifact
    with (
        patch(
            "lamindb.integrations.notion._download_file_to_temp_path",
            return_value="mock-notion-a.pdf",
        ),
        patch(
            "lamindb.integrations.notion.ln.Artifact", return_value=record
        ) as Artifact,
        patch("lamindb.integrations.notion.os.remove"),
    ):
        detail = (
            f"{_short_file_source('https://example.com/a.pdf')} <- row-1:Attachment"
        )
        out = _ensure_artifacts(
            {"https://example.com/a.pdf"},
            transfer_details_by_url={"https://example.com/a.pdf": detail},
            report=report,
        )
    Artifact.assert_called_once_with(
        "mock-notion-a.pdf", key=_artifact_key_from_url("https://example.com/a.pdf")
    )
    assert out["https://example.com/a.pdf"] is created_artifact
    assert report.created_artifacts == [detail]


def test_flatten_created_time(page_props):
    assert _flatten(page_props["Created"]) == "2024-01-10T08:00:00.000Z"


def test_flatten_last_edited_time(page_props):
    assert _flatten(page_props["Edited"]) == "2024-01-15T10:30:00.000Z"


def test_flatten_rollup_returns_none(page_props):
    assert _flatten(page_props["Rollup"]) is None


def test_flatten_formula_returns_none(page_props):
    assert _flatten(page_props["Formula"]) is None


def test_flatten_unknown_type_returns_none():
    assert _flatten({"type": "future_widget", "future_widget": {"value": 99}}) is None


def test_flatten_empty_dict():
    assert _flatten({}) is None


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------


def test_init_raises_on_empty_token():
    with patch("httpx.Client"):
        with pytest.raises(ValueError, match="access token"):
            _NotionReader(token="")  # noqa: S106


def test_init_sets_headers():
    with patch("httpx.Client") as MockSession:
        MockSession.return_value = MagicMock()
        client = _NotionReader(token="tok")  # noqa: S106
    headers = client.s.headers.update.call_args[0][0]
    assert headers["Authorization"] == "Bearer tok"
    assert headers["Notion-Version"] == API_VERSION
    assert headers["Content-Type"] == "application/json"


# ---------------------------------------------------------------------------
# _call — errors and request shape
# ---------------------------------------------------------------------------


def test_call_sends_timeout(reader):
    reader.s.request.return_value = _make_response(DB)
    reader.data_sources("db")
    assert reader.s.request.call_args[1]["timeout"] == 30


def test_401_raises_permission_error(reader):
    reader.s.request.return_value = _make_response({}, 401)
    with pytest.raises(PermissionError, match="token"):
        reader.data_sources("db")


def test_404_mentions_sharing(reader):
    reader.s.request.return_value = _make_response({}, 404)
    with pytest.raises(LookupError, match="not shared with this connection"):
        reader.data_sources("db")


def test_500_falls_through_to_raise_for_status(reader):
    r500 = _make_response({}, 500)
    reader.s.request.return_value = r500
    with patch("lamindb.integrations.notion.time.sleep") as sleep:
        with pytest.raises(RuntimeError, match="gave up after retries"):
            reader.data_sources("db")
    assert [c.args[0] for c in sleep.call_args_list] == [1, 2, 4, 8, 16, 32]
    r500.raise_for_status.assert_not_called()


# ---------------------------------------------------------------------------
# data source resolution
# ---------------------------------------------------------------------------


def test_data_sources_returns_list(reader):
    reader.s.request.return_value = _make_response(DB)
    assert reader.data_sources("db") == [{"id": DS_ID, "name": "Default"}]


def test_resolve_hits_databases_endpoint(reader):
    reader.s.request.return_value = _make_response(DB)
    assert reader._resolve("db-1") == DS_ID
    method, url = reader.s.request.call_args[0]
    assert method == "GET"
    assert url == f"{BASE}/databases/db-1"


def test_resolve_is_cached(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(DS)]
    reader.columns("db-1")
    reader.columns("db-1")
    # resolution and schema are both cached
    assert reader.s.request.call_count == 2


def test_resolve_no_data_sources_raises(reader):
    reader.s.request.return_value = _make_response(
        {"object": "database", "data_sources": []}
    )
    with pytest.raises(LookupError, match="No data sources"):
        reader._resolve("db-1")


def test_resolve_multi_source_warns_and_picks_first(reader):
    multi = {
        "data_sources": [
            {"id": "ds-a", "name": "Alpha"},
            {"id": "ds-b", "name": "Beta"},
        ]
    }
    reader.s.request.return_value = _make_response(multi)
    with patch("lamindb.integrations.notion.logger") as log:
        assert reader._resolve("db-1") == "ds-a"
    msg = log.warning.call_args[0][0]
    assert "2 data sources" in msg
    assert "Alpha" in msg


# ---------------------------------------------------------------------------
# schema & columns
# ---------------------------------------------------------------------------


def test_columns_returns_name_to_type(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(DS)]
    cols = reader.columns("db-1")
    assert cols["Name"] == "title"
    assert cols["Related"] == "relation"
    assert cols["Rollup"] == "rollup"


def test_columns_hits_data_source_endpoint(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(DS)]
    reader.columns("db-1")
    method, url = reader.s.request.call_args_list[1][0]
    assert method == "GET"
    assert url == f"{BASE}/data_sources/{DS_ID}"


def test_schema_exposes_relation_target(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(DS)]
    schema = reader.schema("db-1")
    assert schema["Related"]["target"] == "other-ds-0000"
    assert schema["Name"]["target"] is None


def test_schema_exposes_dual_property(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(DS)]
    schema = reader.schema("db-1")
    assert schema["Related"]["dual"]["synced_property_name"] == "people"
    assert schema["Name"]["dual"] is None  # not a relation


def test_schema_exposes_select_and_multiselect_choices(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(DS)]
    schema = reader.schema("db-1")
    assert schema["Priority"]["choices"] == ["High"]
    assert schema["Tags"]["choices"] == ["python"]


def test_schema_single_property_relation_has_no_dual(reader):
    single = {
        "properties": {
            "Owner": {
                "id": "o",
                "type": "relation",
                "relation": {
                    "data_source_id": "ds-z",
                    "type": "single_property",
                    "single_property": {},
                },
            }
        }
    }
    reader.s.request.side_effect = [_make_response(DB), _make_response(single)]
    schema = reader.schema("db-1")
    assert schema["Owner"]["target"] == "ds-z"
    assert schema["Owner"]["dual"] is None


# ---------------------------------------------------------------------------
# rows — relation values are always UUIDs
# ---------------------------------------------------------------------------


def test_rows_single_page(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    rows = reader.rows("db-1")
    assert len(rows) == 1
    assert rows[0]["notion_id"] == "page-abc-1234"
    assert rows[0]["Name"] == "Row One"
    assert rows[0]["Score"] == 7.5


def test_rows_posts_to_query_endpoint(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    reader.rows("db-1")
    method, url = reader.s.request.call_args_list[1][0]
    assert method == "POST"
    assert url == f"{BASE}/data_sources/{DS_ID}/query"
    assert reader.s.request.call_args_list[1][1]["json"]["page_size"] == 100


def test_rows_pagination_follows_cursor(reader):
    p1 = _make_response(
        {
            "results": [{"id": "p1", "last_edited_time": "t1", "properties": {}}],
            "has_more": True,
            "next_cursor": "cursor-abc",
        }
    )
    p2 = _make_response(
        {
            "results": [{"id": "p2", "last_edited_time": "t2", "properties": {}}],
            "has_more": False,
            "next_cursor": None,
        }
    )
    reader.s.request.side_effect = [_make_response(DB), p1, p2]
    rows = reader.rows("db-1")
    assert [r["notion_id"] for r in rows] == ["p1", "p2"]
    assert reader.s.request.call_args_list[2][1]["json"]["start_cursor"] == "cursor-abc"


def test_rows_empty_returns_empty_list(reader):
    empty = _make_response({"results": [], "has_more": False, "next_cursor": None})
    reader.s.request.side_effect = [_make_response(DB), empty]
    assert reader.rows("db-1") == []


def test_rows_first_two_keys_are_page_level(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    rows = reader.rows("db-1")
    assert list(rows[0])[:3] == ["notion_id", "created_time", "last_edited_time"]


def test_rows_page_fields_win_over_same_named_property(reader):
    clash = {
        "results": [
            {
                "id": "real-page-id",
                "last_edited_time": "2024-05-01T00:00:00Z",
                "properties": {
                    "notion_id": {
                        "type": "rich_text",
                        "rich_text": [{"plain_text": "user-value"}],
                    }
                },
            }
        ],
        "has_more": False,
        "next_cursor": None,
    }
    reader.s.request.side_effect = [_make_response(DB), _make_response(clash)]
    rows = reader.rows("db-1")
    assert rows[0]["notion_id"] == "real-page-id"


def test_rows_preserves_property_names_verbatim(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    rows = reader.rows("db-1")
    for name in ("Name", "Email", "Tags", "Due", "Related", "Author"):
        assert name in rows[0]


def test_rows_never_contain_titles(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    rows = reader.rows("db-1")
    assert rows[0]["Related"] == ["related-page-id-111"]


def test_drop_by_property_name(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    rows = reader.rows("db-1", drop={"Notes", "Related"})
    assert "Notes" not in rows[0]
    assert "Related" not in rows[0]
    assert "Name" in rows[0]


def test_drop_by_property_type(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    rows = reader.rows("db-1", drop={"created_by", "last_edited_by"})
    assert "Author" not in rows[0]
    assert "LastEditor" not in rows[0]


def test_drop_does_not_remove_page_level_keys(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(PAGES)]
    rows = reader.rows("db-1", drop={"notion_id", "last_edited_time"})
    assert rows[0]["notion_id"] == "page-abc-1234"


# ---------------------------------------------------------------------------
# relation_titles — display-only side table
# ---------------------------------------------------------------------------


def test_relation_titles_returns_uuid_to_name(reader):
    reader.s.request.side_effect = [
        _make_response(DB),
        _make_response(DS),
        _make_response(ORG_PAGES),
    ]
    titles = reader.relation_titles("db-1")
    assert titles["related-page-id-111"] == "Deepmind"


def test_relation_titles_empty_when_target_unshared(reader):
    reader.s.request.side_effect = [
        _make_response(DB),
        _make_response(DS),
        _make_response({}, 404),
    ]
    with patch("lamindb.integrations.notion.logger") as log:
        assert reader.relation_titles("db-1") == {}
    assert "not shared" in log.warning.call_args[0][0]


def test_relation_titles_honours_drop(reader):
    reader.s.request.side_effect = [_make_response(DB), _make_response(DS)]
    # Related is the only relation; dropping it means no title queries at all
    assert reader.relation_titles("db-1", drop={"Related"}) == {}
    assert reader.s.request.call_count == 2


def test_title_map_is_cached(reader):
    reader.s.request.side_effect = [_make_response(ORG_PAGES)]
    reader.title_map("ds-x")
    reader.title_map("ds-x")
    assert reader.s.request.call_count == 1


def test_page_title_found_by_type_not_name(reader):
    unnamed = {
        "results": [
            {
                "id": "p-1",
                "properties": {
                    "": {"type": "title", "title": [{"plain_text": "No Name Prop"}]}
                },
            }
        ],
        "has_more": False,
    }
    reader.s.request.side_effect = [_make_response(unnamed)]
    assert reader.title_map("ds-y") == {"p-1": "No Name Prop"}


# ---------------------------------------------------------------------------
# to_json
# ---------------------------------------------------------------------------


def test_to_json_structure(reader):
    reader.s.request.side_effect = [
        _make_response(DB),
        _make_response(DS),
        _make_response(ORG_PAGES),
        _make_response(PAGES),
    ]
    doc = json.loads(reader.to_json("db-1"))
    assert doc["database_id"] == "db-1"
    assert doc["columns"]["Name"] == "title"
    assert doc["titles"]["related-page-id-111"] == "Deepmind"
    assert len(doc["rows"]) == 1
    assert doc["rows"][0]["notion_id"] == "page-abc-1234"
    # rows keep UUIDs, never titles
    assert doc["rows"][0]["Related"] == ["related-page-id-111"]


def test_to_json_writes_file(reader, tmp_path):
    reader.s.request.side_effect = [
        _make_response(DB),
        _make_response(DS),
        _make_response(ORG_PAGES),
        _make_response(PAGES),
    ]
    out = tmp_path / "out.json"
    reader.to_json("db-1", str(out))
    assert json.loads(out.read_text(encoding="utf-8"))["database_id"] == "db-1"


def test_rows_limit_caps_page_size(reader):
    reader.s.request.side_effect = [
        _make_response(DB),
        _make_response(
            {"results": [{"id": "p1", "properties": {}}], "has_more": False}
        ),
    ]
    reader.rows("db-1", limit=5)
    assert reader.s.request.call_args_list[1][1]["json"]["page_size"] == 5


def test_rows_limit_stops_paginating(reader):
    big = _make_response(
        {
            "results": [{"id": f"p{i}", "properties": {}} for i in range(100)],
            "has_more": True,
            "next_cursor": "c",
        }
    )
    reader.s.request.side_effect = [_make_response(DB), big]
    rows = reader.rows("db-1", limit=100)
    assert len(rows) == 100
    assert reader.s.request.call_count == 2  # did not follow the cursor


def test_rows_limit_truncates_overshoot(reader):
    over = _make_response(
        {
            "results": [{"id": f"p{i}", "properties": {}} for i in range(3)],
            "has_more": False,
        }
    )
    reader.s.request.side_effect = [_make_response(DB), over]
    assert len(reader.rows("db-1", limit=2)) == 2


def test_rows_no_limit_paginates_fully(reader):
    p1 = _make_response(
        {
            "results": [{"id": "a", "properties": {}}],
            "has_more": True,
            "next_cursor": "c",
        }
    )
    p2 = _make_response({"results": [{"id": "b", "properties": {}}], "has_more": False})
    reader.s.request.side_effect = [_make_response(DB), p1, p2]
    assert len(reader.rows("db-1")) == 2


def test_page_markdown_exports_nested_blocks(reader):
    reader.s.request.side_effect = [
        _make_response(
            {
                "results": [
                    {
                        "id": "h1",
                        "type": "heading_2",
                        "has_children": False,
                        "heading_2": {"rich_text": [{"plain_text": "Overview"}]},
                    },
                    {
                        "id": "todo-1",
                        "type": "to_do",
                        "has_children": False,
                        "to_do": {
                            "checked": True,
                            "rich_text": [{"plain_text": "Ship integration"}],
                        },
                    },
                    {
                        "id": "list-1",
                        "type": "bulleted_list_item",
                        "has_children": True,
                        "bulleted_list_item": {
                            "rich_text": [{"plain_text": "Milestones"}]
                        },
                    },
                ],
                "has_more": False,
            }
        ),
        _make_response(
            {
                "results": [
                    {
                        "id": "child-1",
                        "type": "paragraph",
                        "has_children": False,
                        "paragraph": {"rich_text": [{"plain_text": "v1 in October"}]},
                    }
                ],
                "has_more": False,
            }
        ),
    ]
    markdown = reader.page_markdown("page-1")
    assert "## Overview" in markdown
    assert "- [x] Ship integration" in markdown
    assert "- Milestones" in markdown
    assert "  v1 in October" in markdown


def test_page_markdown_exports_toggle_as_details_html(reader):
    reader.s.request.side_effect = [
        _make_response(
            {
                "results": [
                    {
                        "id": "toggle-1",
                        "type": "toggle",
                        "has_children": True,
                        "toggle": {"rich_text": [{"plain_text": "Details"}]},
                    }
                ],
                "has_more": False,
            }
        ),
        _make_response(
            {
                "results": [
                    {
                        "id": "paragraph-1",
                        "type": "paragraph",
                        "has_children": False,
                        "paragraph": {"rich_text": [{"plain_text": "Inner note"}]},
                    }
                ],
                "has_more": False,
            }
        ),
    ]
    markdown = reader.page_markdown("page-1")
    assert "<details>" in markdown
    assert "<summary>Details</summary>" in markdown
    assert "<p>" in markdown
    assert "  Inner note" in markdown
    assert "</p>" in markdown
    assert "</details>" in markdown


# ---------------------------------------------------------------------------
# _NotionSyncer
# ---------------------------------------------------------------------------


@pytest.fixture()
def syncer(monkeypatch):
    monkeypatch.setenv("NOTION_TOKEN", "env-token")
    with patch("httpx.Client") as MockSession:
        MockSession.return_value = MagicMock()
        return _NotionSyncer()


def _fake_rec_type(name: str, features: list[str]):
    schema = type(
        "Schema", (), {"members": [type("F", (), {"name": f}) for f in features]}
    )
    return type("RecordType", (), {"name": name, "schema": schema})()


def _ts(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def _fake_record(updated_at: str | None):
    updated = _ts(updated_at) if updated_at else None
    return type("Record", (), {"updated_at": updated, "created_at": updated})()


def test_attach_page_markdown_adds_readme_only_when_changed():
    record = MagicMock()
    record.notes = "old notes"
    saved_block = MagicMock()
    block_factory = MagicMock()
    block_factory.save.return_value = saved_block
    with patch(
        "lamindb.integrations.notion.ln.models.RecordBlock", return_value=block_factory
    ) as RecordBlock:
        _attach_page_markdown(record, "new notes")
    RecordBlock.assert_called_once_with(
        record=record, content="new notes", kind="readme"
    )
    record.ablocks.add.assert_called_once_with(saved_block, bulk=False)

    record.ablocks.add.reset_mock()
    with patch("lamindb.integrations.notion.ln.models.RecordBlock") as RecordBlock:
        _attach_page_markdown(record, "old notes")
    RecordBlock.assert_not_called()
    record.ablocks.add.assert_not_called()


def test_write_attaches_page_markdown_to_records():
    rec = MagicMock()
    rec.notes = None
    rows = [{"notion_id": "page-1", "Name": "A"}]
    by_id = {"page-1": rec}
    reader = MagicMock()
    reader.page_markdown.return_value = "# Notes\n\nhello"
    rec_type = _fake_rec_type("People", ["Name"])
    saved_block = MagicMock()
    block_factory = MagicMock()
    block_factory.save.return_value = saved_block

    with patch(
        "lamindb.integrations.notion.ln.models.RecordBlock", return_value=block_factory
    ):
        stats = _write(
            reader,
            rows,
            rec_type,
            spec={"Name": {"type": "title"}},
            by_id=by_id,
        )

    rec.features.set_values.assert_called_once()
    reader.page_markdown.assert_called_once_with("page-1")
    rec.ablocks.add.assert_called_once_with(saved_block, bulk=False)
    assert stats == {"records": 1, "pending": 0}


def test_syncer_init_raises_without_token(monkeypatch):
    monkeypatch.delenv("NOTION_TOKEN", raising=False)
    with pytest.raises(ValueError, match="NOTION_TOKEN"):
        with patch("httpx.Client") as MockSession:
            MockSession.return_value = MagicMock()
            _NotionSyncer()


def test_import_pages_requires_parents(syncer):
    with pytest.raises(ValueError, match="parents is required"):
        syncer.import_pages([])


def test_resolve_record_type_dry_run_reports_create_record_types(syncer):
    db_id = "3b2d2040-857e-4feb-bb68-d2bec9d6ba09"
    report = SyncReport()
    with (
        patch.object(
            syncer.reader,
            "_call",
            return_value={"title": [{"plain_text": "Website analytics"}]},
        ),
        patch.object(
            syncer.reader, "columns", return_value={"Name": "title", "Score": "number"}
        ),
        patch.object(
            syncer,
            "_database_feature_plan",
            return_value=[("Name", "str", str)],
        ),
        patch.object(syncer, "_plan_or_create_db_metadata"),
        patch("lamindb.integrations.notion.ln.Record") as Record,
    ):
        qs = MagicMock()
        qs.count.return_value = 0
        Record.filter.return_value = qs
        rec_type = syncer._resolve_record_type(db_id, apply=False, report=report)
    assert rec_type is None
    assert report.create_record_types == ["Website analytics"]


def test_resolve_record_type_creates_type_when_missing(syncer):
    db_id = "3b2d2040-857e-4feb-bb68-d2bec9d6ba09"
    report = SyncReport()
    created_type = object()
    with (
        patch.object(
            syncer.reader,
            "_call",
            return_value={"title": [{"plain_text": "Website analytics"}]},
        ),
        patch.object(
            syncer.reader, "columns", return_value={"Name": "title", "Score": "number"}
        ),
        patch.object(
            syncer,
            "_database_feature_plan",
            return_value=[("Name", "str", str)],
        ),
        patch("lamindb.integrations.notion.ln.Record") as Record,
        patch.object(
            syncer, "_create_record_type", return_value=created_type
        ) as create,
    ):
        qs = MagicMock()
        qs.count.return_value = 0
        Record.filter.return_value = qs
        rec_type = syncer._resolve_record_type(db_id, apply=True, report=report)
    assert rec_type is created_type
    create.assert_called_once_with(
        db_id, "Website analytics", None, None, report=report
    )
    assert report.created_record_types == ["Website analytics"]


def test_resolve_record_type_apply_adds_missing_features_to_existing_schema(syncer):
    db_id = "3b2d2040-857e-4feb-bb68-d2bec9d6ba09"
    report = SyncReport()
    existing_feature = type("Feature", (), {"name": "name"})()
    missing_feature = type("Feature", (), {"name": "page_views"})()
    rec_type = MagicMock()
    rec_type.name = "Website analytics"
    rec_type.description = None
    rec_type._aux = None
    rec_type.schema = MagicMock()
    rec_type.schema.members = [existing_feature]

    with (
        patch.object(
            syncer.reader,
            "_call",
            return_value={"title": [{"plain_text": "Website analytics"}]},
        ),
        patch.object(
            syncer.reader,
            "columns",
            return_value={"name": "title", "page_views": "number"},
        ),
        patch.object(
            syncer,
            "_database_feature_plan",
            return_value=[("name", "str", str), ("page_views", "num", "num")],
        ),
        patch.object(
            syncer,
            "_plan_or_create_db_metadata",
            return_value=(
                object(),
                [existing_feature, missing_feature],
                rec_type.schema,
            ),
        ) as plan_or_create,
        patch("lamindb.integrations.notion.ln.Record") as Record,
    ):
        qs = MagicMock()
        qs.count.return_value = 1
        qs.one.return_value = rec_type
        Record.filter.return_value = qs
        resolved = syncer._resolve_record_type(db_id, apply=True, report=report)

    assert resolved is rec_type
    plan_or_create.assert_called_once()
    rec_type.schema.add.assert_called_once_with([missing_feature])


def test_resolve_record_type_dry_run_plans_missing_features_for_existing_schema(syncer):
    db_id = "3b2d2040-857e-4feb-bb68-d2bec9d6ba09"
    report = SyncReport()
    rec_type = MagicMock()
    rec_type.name = "Website analytics"
    rec_type.description = None
    rec_type._aux = None
    rec_type.schema = MagicMock()

    with (
        patch.object(
            syncer.reader,
            "_call",
            return_value={"title": [{"plain_text": "Website analytics"}]},
        ),
        patch.object(
            syncer.reader,
            "columns",
            return_value={"name": "title", "page_views": "number"},
        ),
        patch.object(
            syncer,
            "_database_feature_plan",
            return_value=[("name", "str", str), ("page_views", "num", "num")],
        ),
        patch.object(syncer, "_plan_or_create_db_metadata") as plan_or_create,
        patch("lamindb.integrations.notion.ln.Record") as Record,
    ):
        qs = MagicMock()
        qs.count.return_value = 1
        qs.one.return_value = rec_type
        Record.filter.return_value = qs
        resolved = syncer._resolve_record_type(db_id, apply=False, report=report)

    assert resolved is rec_type
    plan_or_create.assert_called_once()
    assert plan_or_create.call_args.kwargs["apply"] is False


def test_plan_metadata_uses_schema_members_when_feature_type_missing(syncer):
    report = SyncReport(apply=False)
    feature_plan = [
        ("summary", "str", str),
        ("interaction", "list[str]", list[str]),
        ("website", "str", str),
    ]
    schema = MagicMock()
    summary_feature = type("Feature", (), {"name": "summary", "dtype_as_str": "str"})()
    interaction_feature = type(
        "Feature", (), {"name": "interaction", "dtype_as_str": "list[ULabel]"}
    )()
    schema.members.all.return_value = [summary_feature, interaction_feature]
    schema.members.filter.return_value = [summary_feature, interaction_feature]

    with (
        patch("lamindb.integrations.notion.ln.Feature") as Feature,
        patch("lamindb.integrations.notion.ln.Schema") as Schema,
    ):
        feature_type_qs = MagicMock()
        feature_type_qs.count.return_value = 0
        feature_type_qs.one_or_none.return_value = None
        Feature.filter.return_value = feature_type_qs

        schema_qs = MagicMock()
        schema_qs.count.return_value = 1
        schema_qs.one_or_none.return_value = schema
        Schema.filter.return_value = schema_qs

        feature_type, features, returned_schema = syncer._plan_or_create_db_metadata(
            "Organizations",
            feature_plan,
            apply=False,
            report=report,
        )

    assert feature_type is None
    assert returned_schema is schema
    assert {feature.name for feature in features} == {"summary", "interaction"}
    assert report.create_feature_types == ["Organizations"]
    assert report.update_schemas == ["Organizations"]
    assert report.create_features == ["Organizations / website: str"]
    assert report.update_features == [
        "Organizations / summary: str",
        "Organizations / interaction: list[str]",
    ]


def test_plan_metadata_apply_updates_existing_schema_features_to_feature_type(syncer):
    report = SyncReport(apply=True)
    feature_plan = [("summary", "str", str)]
    schema = MagicMock()
    existing_feature = MagicMock()
    existing_feature.name = "summary"
    existing_feature.type_id = None
    schema.members.all.return_value = [existing_feature]
    schema.members.filter.return_value = [existing_feature]

    with (
        patch("lamindb.integrations.notion.ln.Feature") as Feature,
        patch("lamindb.integrations.notion.ln.Schema") as Schema,
    ):
        feature_type = MagicMock()
        feature_type.id = 42
        feature_type_qs = MagicMock()
        feature_type_qs.count.return_value = 0
        feature_type_qs.one_or_none.return_value = None
        Feature.filter.return_value = feature_type_qs
        Feature.return_value.save.return_value = feature_type

        schema_qs = MagicMock()
        schema_qs.count.return_value = 1
        schema_qs.one_or_none.return_value = schema
        Schema.filter.return_value = schema_qs

        syncer._plan_or_create_db_metadata(
            "Organizations",
            feature_plan,
            apply=True,
            report=report,
        )

    Feature.assert_called_once_with(name="Organizations", is_type=True)
    existing_feature.save.assert_called_once_with(update_fields=["type"])
    assert existing_feature.type is feature_type
    assert report.updated_features == ["Organizations / summary: str"]


def test_update_features_report_existing_lamin_dtype_labels(syncer):
    report = SyncReport(apply=False)
    feature_plan = [
        ("business_type", "list[str]", list[str]),
        ("summary", "str", str),
        ("person", "list[str]", list[str]),
        ("interaction", "list[str]", list[str]),
        ("website", "str", str),
    ]
    schema = MagicMock()

    def _feature(name: str, dtype: str):
        feature = MagicMock()
        feature.name = name
        feature.dtype_as_str = dtype
        feature.type_id = None
        return feature

    business_type = _feature("business_type", "list[ULabel]")
    summary = _feature("summary", "str")
    person = _feature("person", "list[People]")
    interaction = _feature("interaction", "list[ULabel]")
    schema.members.all.return_value = [business_type, summary, person, interaction]
    schema.members.filter.return_value = [business_type, summary, person, interaction]

    with (
        patch("lamindb.integrations.notion.ln.Feature") as Feature,
        patch("lamindb.integrations.notion.ln.Schema") as Schema,
    ):
        feature_type_qs = MagicMock()
        feature_type_qs.count.return_value = 0
        feature_type_qs.one_or_none.return_value = None
        Feature.filter.return_value = feature_type_qs

        schema_qs = MagicMock()
        schema_qs.count.return_value = 1
        schema_qs.one_or_none.return_value = schema
        Schema.filter.return_value = schema_qs

        syncer._plan_or_create_db_metadata(
            "Organizations",
            feature_plan,
            apply=False,
            report=report,
        )

    assert report.update_features == [
        "Organizations / business_type: list[str]",
        "Organizations / summary: str",
        "Organizations / person: list[str]",
        "Organizations / interaction: list[str]",
    ]
    assert report.create_features == ["Organizations / website: str"]


def test_create_record_type_uses_title_property_as_schema_index(syncer):
    db_id = "3b2d2040-857e-4feb-bb68-d2bec9d6ba09"
    report = SyncReport()
    columns = {"Display name": "title", "Score": "number"}
    feature_plan = [
        ("Display name", "str", str),
        ("Score", "num", "num"),
    ]
    schema = object()
    with (
        patch.object(syncer.reader, "columns", return_value=columns),
        patch.object(syncer, "_database_feature_plan", return_value=feature_plan),
        patch.object(
            syncer,
            "_plan_or_create_db_metadata",
            return_value=(object(), [], schema),
        ) as plan_or_create,
        patch("lamindb.integrations.notion.ln.Record") as Record,
    ):
        Record.return_value.save.return_value = "created-type"
        rec_type = syncer._create_record_type(
            db_id,
            "Website analytics",
            "Website traffic including page views and unique visitors.",
            "🦆",
            report=report,
        )
    assert rec_type == "created-type"
    plan_or_create.assert_called_once_with(
        "Website analytics",
        feature_plan,
        index_feature_name="Display name",
        record_field_mappings={},
        apply=True,
        report=report,
    )
    Record.assert_called_once_with(
        name="Website analytics",
        description="Website traffic including page views and unique visitors.",
        is_type=True,
        schema=schema,
        _aux={"ei": "🦆"},
    )


def test_record_field_mapping_is_derived_from_internal_notion_property_types(syncer):
    columns = {
        "Created At": "created_time",
        "Edited At": "last_edited_time",
        "Creator": "created_by",
        "Reviewer": "last_edited_by",
        "Name": "title",
    }
    mappings = syncer._record_field_mappings_from_columns(columns)
    assert mappings == {
        "Created At": "created_at",
        "Edited At": "updated_at",
        "Creator": "created_by",
    }


def test_database_emoji_and_description_parsing(syncer):
    payload = {
        "title": [{"plain_text": "Website analytics"}],
        "description": [{"plain_text": "Website traffic including page views."}],
        "icon": {"type": "emoji", "emoji": "🦆"},
    }
    assert syncer._database_title(payload, fallback="db-1") == "Website analytics"
    assert (
        syncer._database_description(payload) == "Website traffic including page views."
    )
    assert syncer._database_emoji(payload) == "🦆"


def test_merge_aux_with_emoji_uses_frontend_convention(syncer):
    assert syncer._merge_aux_with_emoji({"ss": 1}, "🦆") == {"ss": 1, "ei": "🦆"}
    assert syncer._merge_aux_with_emoji({"ss": 1, "ei": "🦆"}, None) == {"ss": 1}


def test_feature_dtype_for_files_maps_to_artifact_list(syncer):
    dtype = syncer._feature_dtype_from_notion_type("files")
    assert syncer._feature_dtype_label_from_notion_type("files") == "list[Artifact]"
    assert getattr(dtype, "__origin__", None) is list
    assert dtype.__args__[0] is ln.Artifact


def test_database_feature_plan_inferrs_multi_select_and_relation_semantics(syncer):
    db_id = "db-1"
    schema_spec = {
        "business_type": {"type": "multi_select", "target": None, "dual": None},
        "person": {
            "type": "relation",
            "target": "ds-people",
            "dual": {"synced_property_name": "people"},
        },
    }
    columns = {"business_type": "multi_select", "person": "relation"}
    people_type = type("PeopleType", (), {"name": "People"})
    people_ulabel_type = type("PeopleULabelType", (), {"name": "People"})

    with (
        patch.object(
            syncer, "_resolve_ulabel_type_by_name_candidates"
        ) as resolve_ulabel,
        patch.object(syncer, "_relation_target_name_candidates") as target_names,
        patch.object(
            syncer, "_resolve_record_type_by_name_candidates"
        ) as resolve_record,
    ):
        resolve_ulabel.return_value = people_ulabel_type
        target_names.return_value = ["People"]
        resolve_record.return_value = people_type
        plan = syncer._database_feature_plan(
            db_id, columns=columns, schema_spec=schema_spec
        )

    assert plan[0][0] == "business_type"
    assert plan[0][1] == "list[People]"
    assert plan[1][0] == "person"
    assert plan[1][1] == "list[People]"
    assert resolve_ulabel.called
    assert resolve_record.called


def test_relation_dtype_with_target_does_not_fallback_to_property_name(syncer):
    with patch.object(
        syncer, "_resolve_record_type_by_name_candidates"
    ) as resolve_record:
        resolve_record.side_effect = [None, None]
        with patch.object(syncer, "_relation_target_name_candidates") as target_names:
            target_names.return_value = []
            dtype_label, dtype = syncer._dtype_from_notion_property(
                "Organizations",
                "software",
                {
                    "type": "relation",
                    "target": "target-ds-id",
                    "dual": {"synced_property_name": "organization"},
                },
            )

    assert dtype_label == "list[str]"
    assert getattr(dtype, "__origin__", None) is list
    assert dtype.__args__[0] is str
    assert resolve_record.call_args_list[0].args[0] == []
    assert "software" in resolve_record.call_args_list[1].args[0]
    assert "Software" in resolve_record.call_args_list[1].args[0]
    assert "organization" not in resolve_record.call_args_list[1].args[0]
    assert "Organizations" not in resolve_record.call_args_list[1].args[0]


def test_relation_dtype_with_target_falls_back_to_property_name_only(syncer):
    reference_type = type("ReferenceType", (), {"name": "References"})
    with patch.object(
        syncer, "_resolve_record_type_by_name_candidates"
    ) as resolve_record:
        resolve_record.side_effect = [None, reference_type]
        with patch.object(syncer, "_relation_target_name_candidates") as target_names:
            target_names.return_value = []
            dtype_label, dtype = syncer._dtype_from_notion_property(
                "Organizations",
                "reference",
                {
                    "type": "relation",
                    "target": "target-ds-id",
                    "dual": {"synced_property_name": "organization"},
                },
            )

    assert dtype_label == "list[References]"
    assert getattr(dtype, "__origin__", None) is list
    assert dtype.__args__[0] is reference_type
    assert resolve_record.call_args_list[0].args[0] == []
    assert "reference" in resolve_record.call_args_list[1].args[0]
    assert "Reference" in resolve_record.call_args_list[1].args[0]
    assert "organization" not in resolve_record.call_args_list[1].args[0]
    assert "Organizations" not in resolve_record.call_args_list[1].args[0]


def test_relation_dtype_with_target_plans_missing_record_type(syncer):
    report = SyncReport(apply=False)
    with patch.object(
        syncer, "_resolve_record_type_by_name_candidates"
    ) as resolve_record:
        resolve_record.side_effect = [None, None]
        with patch.object(syncer, "_relation_target_name_candidates") as target_names:
            target_names.return_value = []
            dtype_label, dtype = syncer._dtype_from_notion_property(
                "Organizations",
                "reference",
                {
                    "type": "relation",
                    "target": "target-ds-id",
                    "dual": {"synced_property_name": "organization"},
                },
                apply=False,
                report=report,
            )

    assert dtype_label == "list[References]"
    assert getattr(dtype, "__origin__", None) is list
    assert dtype.__args__[0] is ln.Record
    assert report.create_record_types == ["References"]


def test_relation_dtype_with_target_prefers_notion_database_name(syncer):
    report = SyncReport(apply=False)
    with patch.object(
        syncer, "_resolve_record_type_by_name_candidates"
    ) as resolve_record:
        resolve_record.side_effect = [None]
        with patch.object(syncer, "_relation_target_name_candidates") as target_names:
            target_names.return_value = [
                "software",
                "Software",
                "softwares",
                "Softwares",
            ]
            dtype_label, dtype = syncer._dtype_from_notion_property(
                "Organizations",
                "software",
                {
                    "type": "relation",
                    "target": "target-ds-id",
                    "dual": {"synced_property_name": "organization"},
                },
                apply=False,
                report=report,
            )

    assert dtype_label == "list[Software]"
    assert getattr(dtype, "__origin__", None) is list
    assert dtype.__args__[0] is ln.Record
    assert report.create_record_types == ["Software"]
    assert "Softwares" not in report.create_record_types


def test_multi_select_ulabel_type_candidates_include_hierarchical_child_name(syncer):
    with patch.object(
        syncer, "_resolve_ulabel_type_by_name_candidates"
    ) as resolve_ulabel:
        resolve_ulabel.return_value = None
        syncer._dtype_from_notion_property(
            "Organizations",
            "interaction",
            {"type": "multi_select", "target": None, "dual": None},
        )
    child_candidates = resolve_ulabel.call_args_list[1].args[0]
    assert "interactions" in child_candidates


def test_multi_select_plans_ulabel_type_and_labels_in_dry_run(syncer):
    report = SyncReport(apply=False)
    with patch.object(
        syncer, "_resolve_ulabel_type_by_name_candidates"
    ) as resolve_ulabel:
        resolve_ulabel.return_value = None
        dtype_label, dtype = syncer._dtype_from_notion_property(
            "Organizations",
            "interaction",
            {
                "type": "multi_select",
                "choices": ["Email", "Call"],
                "target": None,
                "dual": None,
            },
            apply=False,
            report=report,
        )
    assert dtype_label == "list[interactions]"
    assert getattr(dtype, "__origin__", None) is list
    assert dtype.__args__[0] is ln.ULabel
    assert report.create_ulabel_types == [
        "Organizations",
        "Organizations / interactions",
    ]
    assert report.create_ulabels == [
        "Organizations / interactions / Call",
        "Organizations / interactions / Email",
    ]


def test_multi_select_creates_ulabel_type_and_labels_in_apply(syncer):
    report = SyncReport(apply=True)
    parent_type = MagicMock()
    parent_type.name = "Organizations"
    parent_type.id = 41
    created_type = MagicMock()
    created_type.name = "interactions"
    created_type.id = 77
    with (
        patch.object(
            syncer,
            "_resolve_ulabel_type_by_name_candidates",
            side_effect=[parent_type, None],
        ),
        patch("lamindb.integrations.notion.ln.ULabel") as ULabel,
    ):
        ULabel.return_value.save.side_effect = [created_type, MagicMock()]
        ULabel.filter.return_value.values_list.return_value = ["Call"]
        dtype_label, _ = syncer._dtype_from_notion_property(
            "Organizations",
            "interaction",
            {
                "type": "multi_select",
                "choices": ["Email", "Call"],
                "target": None,
                "dual": None,
            },
            apply=True,
            report=report,
        )
    assert dtype_label == "list[interactions]"
    assert report.created_ulabel_types == ["Organizations / interactions"]
    assert report.created_ulabels == ["Organizations / interactions / Email"]


def test_collect_database_ids_falls_back_to_page_on_database_400(syncer):
    page_id = "7283894209c44522a7c79620795d0409"
    request = httpx.Request("GET", f"{BASE}/databases/{page_id}")
    response = httpx.Response(400, request=request)
    db_400 = _make_response({}, 400)
    db_400.raise_for_status.side_effect = httpx.HTTPStatusError(
        "bad request",
        request=request,
        response=response,
    )
    page_ok = _make_response({"object": "page", "id": page_id})
    children = _make_response(
        {
            "results": [
                {"id": "db-1", "type": "child_database", "has_children": False}
            ],
            "has_more": False,
            "next_cursor": None,
        }
    )
    syncer.reader.s.request.side_effect = [db_400, page_ok, children]

    db_ids, parent_pages = syncer._collect_database_ids([page_id])

    assert db_ids == {"db-1"}
    assert parent_pages == {page_id: page_id}
    assert syncer._parent_page_emojis == {page_id: None}


def test_collect_database_ids_stores_parent_page_emoji(syncer):
    page_id = "7283894209c44522a7c79620795d0409"
    request = httpx.Request("GET", f"{BASE}/databases/{page_id}")
    response = httpx.Response(400, request=request)
    db_400 = _make_response({}, 400)
    db_400.raise_for_status.side_effect = httpx.HTTPStatusError(
        "bad request",
        request=request,
        response=response,
    )
    page_ok = _make_response(
        {"object": "page", "id": page_id, "icon": {"type": "emoji", "emoji": "📊"}}
    )
    children = _make_response(
        {
            "results": [
                {"id": "db-1", "type": "child_database", "has_children": False}
            ],
            "has_more": False,
            "next_cursor": None,
        }
    )
    syncer.reader.s.request.side_effect = [db_400, page_ok, children]

    syncer._collect_database_ids([page_id])

    assert syncer._parent_page_emojis == {page_id: "📊"}


def test_collect_database_ids_limit_caps_discovered_children_for_page_parent(syncer):
    page_id = "7283894209c44522a7c79620795d0409"
    request = httpx.Request("GET", f"{BASE}/databases/{page_id}")
    response = httpx.Response(400, request=request)
    db_400 = _make_response({}, 400)
    db_400.raise_for_status.side_effect = httpx.HTTPStatusError(
        "bad request",
        request=request,
        response=response,
    )
    page_ok = _make_response({"object": "page", "id": page_id})
    children = _make_response(
        {
            "results": [
                {"id": "db-1", "type": "child_database", "has_children": False},
                {"id": "db-2", "type": "child_database", "has_children": False},
            ],
            "has_more": False,
            "next_cursor": None,
        }
    )
    syncer.reader.s.request.side_effect = [db_400, page_ok, children]

    db_ids, _ = syncer._collect_database_ids([page_id], limit=1)

    assert db_ids == {"db-1"}


def test_collect_database_ids_limit_zero_skips_child_database_traversal(syncer):
    page_id = "7283894209c44522a7c79620795d0409"
    request = httpx.Request("GET", f"{BASE}/databases/{page_id}")
    response = httpx.Response(400, request=request)
    db_400 = _make_response({}, 400)
    db_400.raise_for_status.side_effect = httpx.HTTPStatusError(
        "bad request",
        request=request,
        response=response,
    )
    page_ok = _make_response({"object": "page", "id": page_id})
    syncer.reader.s.request.side_effect = [db_400, page_ok]

    db_ids, parent_pages = syncer._collect_database_ids([page_id], limit=0)

    assert db_ids == set()
    assert parent_pages == {page_id: page_id}
    assert syncer.reader.s.request.call_count == 2


def test_schema_validation_accepts_exact_property_parity(syncer):
    rec_type = _fake_rec_type("People", ["Name"])
    with patch.object(syncer.reader, "columns", return_value={"Name": "title"}):
        syncer._validate_schema("db-1", rec_type)


def test_schema_validation_checks_full_property_parity(syncer):
    rec_type = _fake_rec_type("People", ["Name", "Extra"])
    with patch.object(
        syncer.reader, "columns", return_value={"Name": "title", "Email": "email"}
    ):
        with pytest.raises(ValueError, match="missing in Lamin schema"):
            syncer._validate_schema("db-1", rec_type)


def test_schema_validation_dry_run_allows_missing_without_extra(syncer):
    rec_type = _fake_rec_type("People", ["Name"])
    with patch.object(
        syncer.reader, "columns", return_value={"Name": "title", "Email": "email"}
    ):
        syncer._validate_schema("db-1", rec_type, apply=False)


def test_resolve_record_type_apply_assigns_parent_type(syncer):
    db_id = "3b2d2040-857e-4feb-bb68-d2bec9d6ba09"
    report = SyncReport()
    rec_type = MagicMock()
    rec_type.name = "Website analytics"
    rec_type.description = None
    rec_type._aux = {"ei": "📊"}
    rec_type.schema = MagicMock()
    rec_type.schema.members = []
    rec_type.type_id = None
    parent_type = type("ParentType", (), {"id": 7})()

    with (
        patch.object(syncer.reader, "columns", return_value={"name": "title"}),
        patch.object(
            syncer,
            "_database_feature_plan",
            return_value=[("name", "str", str)],
        ),
        patch.object(
            syncer,
            "_plan_or_create_db_metadata",
            return_value=(object(), [], rec_type.schema),
        ),
        patch("lamindb.integrations.notion.ln.Record") as Record,
    ):
        qs = MagicMock()
        qs.count.return_value = 1
        qs.one.return_value = rec_type
        Record.filter.return_value = qs
        syncer._resolve_record_type(
            db_id,
            apply=True,
            report=report,
            payload={"title": [{"plain_text": "Website analytics"}]},
            parent_type=parent_type,
        )

    assert rec_type.type is parent_type
    rec_type.save.assert_called_once()


def test_import_pages_dry_run_does_not_write(syncer):
    rec_type = _fake_rec_type("People", ["Name"])
    rows = [
        {"notion_id": "a", "last_edited_time": "2024-01-01T00:00:00Z", "Name": "A"},
        {"notion_id": "b", "last_edited_time": "2024-01-02T00:00:00Z", "Name": "B"},
    ]
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=({"db-1"}, {"parent-id": "Parent"}),
        ) as collect_ids,
        patch.object(syncer, "_resolve_record_type", return_value=rec_type),
        patch.object(syncer, "_validate_schema"),
        patch.object(syncer.reader, "rows", return_value=rows),
        patch.object(syncer.reader, "schema", return_value={}),
        patch("lamindb.integrations.notion._existing_by_ref", return_value={}),
        patch("lamindb.integrations.notion._upsert_all") as upsert_all,
        patch("lamindb.integrations.notion._write") as write,
    ):
        report = syncer.import_pages("parent", apply=False)
    collect_ids.assert_called_once_with(["parent"], limit=None)
    assert report.apply is False
    assert report.message == "Dry run report -- nothing got created"
    assert report.discovered_pages == 4
    assert report.discovered == 2
    assert report.created == 2
    assert report.updated == 0
    assert report.unchanged == 0
    upsert_all.assert_not_called()
    write.assert_not_called()


def test_import_pages_dry_run_counts_rows_for_missing_record_type(syncer):
    rows = [
        {"notion_id": "a", "last_edited_time": "2024-01-01T00:00:00Z", "Name": "A"},
        {"notion_id": "b", "last_edited_time": "2024-01-02T00:00:00Z", "Name": "B"},
    ]
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=({"db-1"}, {"parent-id": "Parent"}),
        ),
        patch.object(syncer, "_resolve_record_type", return_value=None),
        patch.object(syncer.reader, "rows", return_value=rows),
        patch.object(syncer.reader, "schema", return_value={}),
        patch("lamindb.integrations.notion._upsert_all") as upsert_all,
        patch("lamindb.integrations.notion._write") as write,
    ):
        report = syncer.import_pages("parent", apply=False)
    assert report.apply is False
    assert report.message == "Dry run report -- nothing got created"
    assert report.discovered == 2
    assert report.created == 2
    upsert_all.assert_not_called()
    write.assert_not_called()


def test_import_pages_dry_run_reports_pending_file_transfers(syncer):
    rec_type = _fake_rec_type("People", ["Name", "Attachment"])
    rows = [
        {
            "notion_id": "a",
            "last_edited_time": "2024-01-01T00:00:00Z",
            "Name": "A",
            "Attachment": ["https://example.com/a.pdf"],
        }
    ]
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=({"db-1"}, {"parent-id": "Parent"}),
        ),
        patch.object(syncer, "_resolve_record_type", return_value=rec_type),
        patch.object(syncer, "_validate_schema"),
        patch.object(syncer.reader, "rows", return_value=rows),
        patch.object(
            syncer.reader, "schema", return_value={"Attachment": {"type": "files"}}
        ),
        patch("lamindb.integrations.notion._existing_by_ref", return_value={}),
        patch("lamindb.integrations.notion._upsert_all") as upsert_all,
        patch("lamindb.integrations.notion._write") as write,
    ):
        report = syncer.import_pages("parent", apply=False)
    assert report.create_artifacts == [
        f"{_short_file_source('https://example.com/a.pdf')} <- a:Attachment"
    ]
    upsert_all.assert_not_called()
    write.assert_not_called()


def test_import_pages_dry_run_includes_parent_page_type(syncer):
    rows = [{"notion_id": "a", "last_edited_time": "2024-01-01T00:00:00Z", "Name": "A"}]
    rec_type = _fake_rec_type("Website analytics", ["Name"])
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=(
                {"db-1"},
                {"7283894209c44522a7c79620795d0409": "Import metrics"},
            ),
        ),
        patch.object(syncer, "_resolve_record_type", return_value=rec_type),
        patch.object(syncer, "_validate_schema"),
        patch.object(syncer.reader, "rows", return_value=rows),
        patch.object(syncer.reader, "schema", return_value={}),
        patch("lamindb.integrations.notion._existing_by_ref", return_value={}),
        patch("lamindb.integrations.notion._upsert_all") as upsert_all,
        patch("lamindb.integrations.notion._write") as write,
        patch("lamindb.integrations.notion.ln.Record") as Record,
    ):
        qs = MagicMock()
        qs.count.return_value = 0
        Record.filter.return_value = qs
        report = syncer.import_pages("parent", apply=False)
    assert "Import metrics" in report.create_record_types
    upsert_all.assert_not_called()
    write.assert_not_called()


def test_import_pages_report_compacts_database_ids(syncer):
    rec_type = _fake_rec_type("People", ["Name"])
    rows = [{"notion_id": "a", "last_edited_time": "2024-01-01T00:00:00Z", "Name": "A"}]
    db_id = "3b2d2040-857e-4feb-bb68-d2bec9d6ba09"
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=({db_id}, {"parent-id": "Parent"}),
        ),
        patch.object(syncer, "_resolve_record_type", return_value=rec_type),
        patch.object(syncer, "_validate_schema"),
        patch.object(syncer.reader, "rows", return_value=rows),
        patch.object(syncer.reader, "schema", return_value={}),
        patch("lamindb.integrations.notion._existing_by_ref", return_value={}),
        patch("lamindb.integrations.notion._upsert_all", return_value={}),
        patch(
            "lamindb.integrations.notion._write",
            return_value={"records": 0, "pending": 0},
        ),
    ):
        report = syncer.import_pages("parent", apply=False)
    assert report.databases == ["3b2d2040857e4febbb68d2bec9d6ba09"]


def test_import_pages_writes_only_created_or_changed(syncer):
    rec_type = _fake_rec_type("People", ["Name"])
    rows = [
        {"notion_id": "a", "last_edited_time": "2024-01-01T00:00:00Z", "Name": "A"},
        {"notion_id": "b", "last_edited_time": "2024-01-02T00:00:00Z", "Name": "B"},
        {"notion_id": "c", "last_edited_time": "2024-01-03T00:00:00Z", "Name": "C"},
    ]
    existing = {
        "a": _fake_record("2024-01-01T00:00:00Z"),
        "b": _fake_record("2024-01-01T00:00:00Z"),
    }
    after = {
        "a": _fake_record("2024-01-01T00:00:00Z"),
        "b": _fake_record("2024-01-01T00:00:00Z"),
        "c": _fake_record(None),
    }
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=({"db-1"}, {"parent-id": "Parent"}),
        ),
        patch.object(syncer, "_resolve_record_type", return_value=rec_type),
        patch.object(syncer, "_validate_schema"),
        patch.object(syncer.reader, "rows", return_value=rows),
        patch.object(syncer.reader, "schema", return_value={}),
        patch("lamindb.integrations.notion._existing_by_ref", return_value=existing),
        patch("lamindb.integrations.notion._upsert_all", return_value=after),
        patch(
            "lamindb.integrations.notion._write",
            return_value={"records": 2, "pending": 1},
        ) as write,
    ):
        report = syncer.import_pages(["parent"], apply=True)
    write_rows = write.call_args[0][1]
    assert [r["notion_id"] for r in write_rows] == ["b", "c"]
    assert report.created == 1
    assert report.updated == 1
    assert report.unchanged == 1
    assert report.pending_relations == 1


def test_import_pages_passes_limit_to_database_discovery(syncer):
    rec_type = _fake_rec_type("People", ["Name"])
    rows = [{"notion_id": "a", "last_edited_time": "2024-01-01T00:00:00Z", "Name": "A"}]
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=({"db-1"}, {"parent-id": "Parent"}),
        ) as collect_ids,
        patch.object(syncer, "_resolve_record_type", return_value=rec_type),
        patch.object(syncer, "_validate_schema"),
        patch.object(syncer.reader, "rows", return_value=rows),
        patch.object(syncer.reader, "schema", return_value={}),
        patch("lamindb.integrations.notion._existing_by_ref", return_value={}),
        patch("lamindb.integrations.notion._upsert_all", return_value={}),
        patch(
            "lamindb.integrations.notion._write",
            return_value={"records": 0, "pending": 0},
        ),
    ):
        syncer.import_pages("parent", apply=True, limit=1)
    collect_ids.assert_called_once_with(["parent"], limit=1)


def test_import_pages_limit_zero_ingests_only_parent_pages(syncer):
    with (
        patch.object(
            syncer,
            "_collect_database_ids",
            return_value=(set(), {"parent-id": "Parent"}),
        ) as collect_ids,
        patch("lamindb.integrations.notion.ln.Record") as Record,
    ):
        qs = MagicMock()
        qs.count.return_value = 0
        Record.filter.return_value = qs
        report = syncer.import_pages("parent", apply=False, limit=0)
    collect_ids.assert_called_once_with(["parent"], limit=0)
    assert report.discovered_pages == 1
    assert report.databases == []
    assert report.discovered == 0
    assert report.create_record_types == ["Parent"]


def test_upsert_all_populates_created_and_updated_from_notion():
    rows = [
        {
            "notion_id": "page-1",
            "name": "A",
            "created_time": "2024-01-01T08:00:00Z",
            "last_edited_time": "2024-01-03T10:00:00Z",
        }
    ]
    rec = MagicMock()
    rec.save.return_value = rec
    with (
        patch("lamindb.integrations.notion._existing_by_ref", return_value={}),
        patch("lamindb.integrations.notion.ln.Record", return_value=rec),
    ):
        out = _upsert_all(rec_type=object(), rows=rows)
    assert out["page-1"] is rec
    assert rec.created_at == _ts("2024-01-01T08:00:00Z")
    assert rec.updated_at == _ts("2024-01-03T10:00:00Z")


def test_sync_from_notion_delegates_to_syncer_and_prints():
    sync_report = SyncReport(created=1, apply=False)
    with (
        patch("lamindb.integrations.notion._NotionSyncer") as Syncer,
        patch("lamindb.integrations.notion.RICH_CONSOLE.print") as rich_print,
    ):
        Syncer.return_value.import_pages.return_value = sync_report
        report = sync_from_notion(parents=["p1", "p2"], apply=False, limit=3)
    Syncer.assert_called_once_with(token=None)
    Syncer.return_value.import_pages.assert_called_once_with(
        parents=["p1", "p2"], apply=False, limit=3
    )
    rich_print.assert_called_once()
    assert "Dry run: nothing got created." in rich_print.call_args[0][0]
    assert report is sync_report


def test_sync_report_pretty_text_groups_and_labels_metrics():
    report = SyncReport(
        apply=False,
        discovered_pages=7,
        databases=["3b2d2040857e4febbb68d2bec9d6ba09"],
        discovered=5,
        created=5,
        updated=1,
        unchanged=2,
        pending_relations=0,
        failed=0,
        create_record_types=["Website analytics"],
        create_feature_types=["Website analytics"],
        create_schemas=["Website analytics"],
        update_schemas=["Website analytics"],
        create_features=[
            "Website analytics / Name: str",
            "Website analytics / Score: num",
        ],
        update_features=["Website analytics / Existing: str"],
        create_artifacts=[
            f"{_short_file_source('https://example.com/a.pdf')} <- row-1:Attachment",
            f"{_short_file_source('https://example.com/b.pdf')} <- row-1:Attachment",
        ],
    )
    text = report.to_pretty_text()
    assert "Discovered 7 Notion pages." in text
    assert "Dry run: nothing got created." in text
    assert "pass apply=True or --apply on the CLI" in text
    assert "Scope" in text
    assert "Actions" in text
    assert "discovered_databases" not in text
    assert "database_ids" not in text
    assert "discovered_records" not in text
    assert "create_feature_types" in text
    assert "create_schemas" in text
    assert "update_schemas" in text
    assert "update_features" in text
    assert "create_features" in text
    assert "Website analytics / Score: num" in text
    assert "create_artifacts" in text
    assert (
        f"{_short_file_source('https://example.com/a.pdf')} <- row-1:Attachment" in text
    )
    assert "create_records" in text
    assert text.index("create_record_types") < text.index("create_records")


def test_sync_report_rich_render_preserves_bracketed_dtypes():
    report = SyncReport(
        apply=False,
        update_features=[
            "Organizations / interaction: list[str]",
            "Organizations / file: list[Artifact]",
        ],
    )
    console = Console(record=True, force_terminal=False, color_system=None, width=200)
    console.print(report.to_pretty_text(), markup=True, highlight=False)
    rendered = console.export_text()
    assert "list[str]" in rendered
    assert "list[Artifact]" in rendered


def test_sync_from_notion_live_smoke_with_env_token():
    token = os.getenv("NOTION_TOKEN")
    run_live = os.getenv("CI") or os.getenv("LAMINDB_RUN_NOTION_LIVE_TESTS") in {
        "1",
        "true",
        "True",
    }
    if not token or not run_live:
        pytest.skip(
            "Set NOTION_TOKEN and run in CI, or set "
            "LAMINDB_RUN_NOTION_LIVE_TESTS=true for a local live smoke test."
        )
    report = sync_from_notion(
        parents="7283894209c44522a7c79620795d0409",
        token=token,
        apply=False,
    )
    assert report.apply is False
    assert report.discovered_pages >= 1
