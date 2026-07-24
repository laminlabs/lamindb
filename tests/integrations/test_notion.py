"""Unit tests for lamindb.integrations.notion.

No live network calls — all HTTP is mocked at the requests.Session level.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from lamindb.integrations.notion import API_VERSION, BASE, Reader, _flatten

FIXTURES = Path(__file__).parent / "fixtures"


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
    """Reader with requests.Session replaced by a MagicMock."""
    with patch("requests.Session") as MockSession:
        MockSession.return_value = MagicMock()
        client = Reader(token="secret-test-token")  # noqa: S106
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
    with patch("requests.Session"):
        with pytest.raises(ValueError, match="access token"):
            Reader(token="")  # noqa: S106


def test_init_sets_headers():
    with patch("requests.Session") as MockSession:
        MockSession.return_value = MagicMock()
        client = Reader(token="tok")  # noqa: S106
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
    r500.raise_for_status.side_effect = RuntimeError("500 Server Error")
    reader.s.request.return_value = r500
    with pytest.raises(RuntimeError, match="500"):
        reader.data_sources("db")
    r500.raise_for_status.assert_called_once()


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
    assert list(rows[0])[:2] == ["notion_id", "last_edited_time"]


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
