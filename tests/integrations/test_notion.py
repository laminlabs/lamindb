"""Unit tests for lamindb.integrations._notion.

No live network calls — all HTTP is mocked at the requests.Session level.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from lamindb.integrations._notion import Notion, _flatten_property

FIXTURES = Path(__file__).parent / "fixtures"


def _load_fixture(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text())


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_response(
    data: dict, status_code: int = 200, headers: dict | None = None
) -> MagicMock:
    resp = MagicMock()
    resp.status_code = status_code
    resp.ok = status_code < 400
    resp.json.return_value = data
    resp.headers = headers or {}
    return resp


@pytest.fixture()
def notion():
    """Notion client with requests.Session replaced by a MagicMock."""
    with patch("requests.Session") as MockSession:
        mock_session = MagicMock()
        MockSession.return_value = mock_session
        client = Notion(token="secret-test-token")  # noqa: S106
    return client


@pytest.fixture(scope="module")
def page_props():
    return _load_fixture("notion_page.json")["results"][0]["properties"]


@pytest.fixture(scope="module")
def db_fixture():
    return _load_fixture("notion_database.json")


# ---------------------------------------------------------------------------
# _flatten_property — one test per type in the dispatch table
# ---------------------------------------------------------------------------


def test_flatten_title(page_props):
    assert _flatten_property(page_props["Name"]) == "Row One"


def test_flatten_title_empty():
    assert _flatten_property({"type": "title", "title": []}) == ""


def test_flatten_rich_text(page_props):
    assert _flatten_property(page_props["Notes"]) == "some notes"


def test_flatten_rich_text_multipart():
    prop = {
        "type": "rich_text",
        "rich_text": [{"plain_text": "foo"}, {"plain_text": " bar"}],
    }
    assert _flatten_property(prop) == "foo bar"


def test_flatten_email(page_props):
    assert _flatten_property(page_props["Email"]) == "alice@example.com"


def test_flatten_email_null():
    assert _flatten_property({"type": "email", "email": None}) is None


def test_flatten_phone_number(page_props):
    assert _flatten_property(page_props["Phone"]) == "+1-555-0100"


def test_flatten_url(page_props):
    assert _flatten_property(page_props["Website"]) == "https://example.com"


def test_flatten_number(page_props):
    assert _flatten_property(page_props["Score"]) == 7.5


def test_flatten_checkbox_true(page_props):
    assert _flatten_property(page_props["Active"]) is True


def test_flatten_checkbox_false():
    assert _flatten_property({"type": "checkbox", "checkbox": False}) is False


def test_flatten_select(page_props):
    assert _flatten_property(page_props["Priority"]) == "High"


def test_flatten_select_null():
    assert _flatten_property({"type": "select", "select": None}) is None


def test_flatten_status(page_props):
    assert _flatten_property(page_props["State"]) == "Done"


def test_flatten_status_null():
    assert _flatten_property({"type": "status", "status": None}) is None


def test_flatten_multi_select(page_props):
    assert _flatten_property(page_props["Tags"]) == ["python", "data"]


def test_flatten_multi_select_empty():
    assert _flatten_property({"type": "multi_select", "multi_select": []}) == []


def test_flatten_date_returns_start(page_props):
    assert _flatten_property(page_props["Due"]) == "2024-02-01"


def test_flatten_date_null():
    assert _flatten_property({"type": "date", "date": None}) is None


def test_flatten_people(page_props):
    assert _flatten_property(page_props["Owner"]) == ["user-id-abc"]


def test_flatten_relation(page_props):
    assert _flatten_property(page_props["Related"]) == ["related-page-id-111"]


def test_flatten_files(page_props):
    urls = _flatten_property(page_props["Attachment"])
    assert urls == ["https://example.com/doc.pdf", "https://s3.amazonaws.com/img.png"]


def test_flatten_created_time(page_props):
    assert _flatten_property(page_props["Created"]) == "2024-01-10T08:00:00.000Z"


def test_flatten_last_edited_time(page_props):
    assert _flatten_property(page_props["Edited"]) == "2024-01-15T10:30:00.000Z"


def test_flatten_rollup_returns_none(page_props):
    assert _flatten_property(page_props["Rollup"]) is None


def test_flatten_formula_returns_none(page_props):
    assert _flatten_property(page_props["Formula"]) is None


def test_flatten_unknown_returns_none_and_logs_debug(caplog):
    prop = {"type": "future_widget", "future_widget": {"value": 99}}
    with caplog.at_level(logging.DEBUG, logger="lamindb.integrations._notion"):
        result = _flatten_property(prop)
    assert result is None
    assert "future_widget" in caplog.text


# ---------------------------------------------------------------------------
# Notion.__init__ — token resolution
# ---------------------------------------------------------------------------


def test_init_raises_on_empty_token():
    with patch("requests.Session"):
        with pytest.raises(ValueError, match="integration token"):
            Notion(token="")  # noqa: S106


def test_init_sets_auth_and_version_headers():
    with patch("requests.Session") as MockSession:
        mock_session = MagicMock()
        MockSession.return_value = mock_session
        client = Notion(token="tok")  # noqa: S106
    headers = client._session.headers.update.call_args[0][0]
    assert headers["Authorization"] == "Bearer tok"
    assert headers["Notion-Version"] == "2022-06-28"
    assert headers["Content-Type"] == "application/json"


# ---------------------------------------------------------------------------
# get_schema
# ---------------------------------------------------------------------------


def test_get_schema_basic(notion, db_fixture):
    notion._session.request.return_value = _make_response(db_fixture)
    schema = notion.get_schema("b86daf14-2a54-4728-bda2-496c5760d863")

    assert schema["Name"]["id"] == "title"  # stable across renames
    assert schema["Name"]["type"] == "title"
    assert schema["Name"]["relation_database_id"] is None
    assert schema["Name"]["options"] is None


def test_get_schema_relation_has_database_id(notion, db_fixture):
    notion._session.request.return_value = _make_response(db_fixture)
    schema = notion.get_schema("b86daf14-2a54-4728-bda2-496c5760d863")

    assert schema["Related"]["type"] == "relation"
    assert schema["Related"]["relation_database_id"] == "other-db-id-0000"


def test_get_schema_select_has_options(notion, db_fixture):
    notion._session.request.return_value = _make_response(db_fixture)
    schema = notion.get_schema("b86daf14-2a54-4728-bda2-496c5760d863")

    assert schema["Priority"]["options"] == ["High", "Medium", "Low"]
    assert schema["State"]["options"] == ["Todo", "In Progress", "Done"]
    assert schema["Tags"]["options"] == ["python", "data"]


def test_get_schema_404_mentions_sharing(notion):
    notion._session.request.return_value = _make_response({}, status_code=404)
    with pytest.raises(LookupError, match="not be shared with the integration"):
        notion.get_schema("missing-db")


# ---------------------------------------------------------------------------
# query_database — pagination
# ---------------------------------------------------------------------------


def test_query_database_single_page(notion):
    data = _load_fixture("notion_page.json")
    notion._session.request.return_value = _make_response(data)
    pages = notion.query_database("b86daf14-2a54-4728-bda2-496c5760d863")
    assert len(pages) == 1
    assert pages[0]["id"] == "page-abc-1234"


def test_query_database_pagination(notion):
    page1 = _make_response(
        {
            "results": [
                {
                    "id": "p1",
                    "last_edited_time": "2024-01-01T00:00:00Z",
                    "properties": {},
                }
            ],
            "has_more": True,
            "next_cursor": "cursor-abc",
        }
    )
    page2 = _make_response(
        {
            "results": [
                {
                    "id": "p2",
                    "last_edited_time": "2024-01-02T00:00:00Z",
                    "properties": {},
                }
            ],
            "has_more": False,
            "next_cursor": None,
        }
    )
    notion._session.request.side_effect = [page1, page2]

    pages = notion.query_database("some-db")

    assert len(pages) == 2
    assert pages[0]["id"] == "p1"
    assert pages[1]["id"] == "p2"
    # second call must carry the cursor
    second_call_body = notion._session.request.call_args_list[1][1]["json"]
    assert second_call_body["start_cursor"] == "cursor-abc"


def test_query_database_since_adds_filter(notion):
    notion._session.request.return_value = _make_response(
        {"results": [], "has_more": False, "next_cursor": None}
    )
    notion.query_database("db-id", since="2024-06-01T00:00:00Z")
    body = notion._session.request.call_args[1]["json"]
    assert body["filter"]["last_edited_time"]["on_or_after"] == "2024-06-01T00:00:00Z"


def test_query_database_404_mentions_sharing(notion):
    notion._session.request.return_value = _make_response({}, status_code=404)
    with pytest.raises(LookupError, match="not be shared with the integration"):
        notion.query_database("missing-db")


def test_query_database_401_mentions_token(notion):
    notion._session.request.return_value = _make_response({}, status_code=401)
    with pytest.raises(PermissionError, match="token"):
        notion.query_database("some-db")


# ---------------------------------------------------------------------------
# to_dataframe
# ---------------------------------------------------------------------------


def test_to_dataframe_columns(notion):
    data = _load_fixture("notion_page.json")
    notion._session.request.return_value = _make_response(data)
    df = notion.to_dataframe("b86daf14-2a54-4728-bda2-496c5760d863")

    assert list(df.columns[:2]) == ["notion_id", "last_edited_time"]
    assert df.loc[0, "notion_id"] == "page-abc-1234"
    assert df.loc[0, "Name"] == "Row One"
    assert df.loc[0, "Score"] == 7.5


def test_to_dataframe_preserves_property_names(notion):
    data = _load_fixture("notion_page.json")
    notion._session.request.return_value = _make_response(data)
    df = notion.to_dataframe("db-id")
    # user property names must not be renamed
    for col in ("Name", "Email", "Tags", "Due", "Related"):
        assert col in df.columns


def test_to_dataframe_since_forwarded(notion):
    # Verify that since= is threaded through to query_database's HTTP body.
    db_fixture = _load_fixture("notion_database.json")
    # second call is get_schema for empty-result column inference
    notion._session.request.side_effect = [
        _make_response({"results": [], "has_more": False, "next_cursor": None}),
        _make_response(db_fixture),
    ]
    notion.to_dataframe("db-id", since="2024-06-01T00:00:00Z")
    first_call_body = notion._session.request.call_args_list[0][1]["json"]
    assert (
        first_call_body["filter"]["last_edited_time"]["on_or_after"]
        == "2024-06-01T00:00:00Z"
    )


def test_to_dataframe_empty_returns_schema_columns(notion):
    db_fixture = _load_fixture("notion_database.json")
    notion._session.request.side_effect = [
        _make_response({"results": [], "has_more": False, "next_cursor": None}),
        _make_response(db_fixture),
    ]
    df = notion.to_dataframe("db-id")
    assert len(df) == 0
    assert "notion_id" in df.columns
    assert "last_edited_time" in df.columns
    assert "Name" in df.columns


def test_to_dataframe_page_fields_win_over_same_named_property(notion):
    # A database that has a user property literally named "notion_id" should
    # not clobber the page-level id in the resulting DataFrame.
    page_data = {
        "results": [
            {
                "object": "page",
                "id": "real-page-id",
                "last_edited_time": "2024-05-01T00:00:00Z",
                "properties": {
                    "notion_id": {
                        "type": "rich_text",
                        "rich_text": [{"plain_text": "user-supplied-value"}],
                    }
                },
            }
        ],
        "has_more": False,
        "next_cursor": None,
    }
    notion._session.request.return_value = _make_response(page_data)
    df = notion.to_dataframe("db-id")
    # page-level field must win
    assert df.loc[0, "notion_id"] == "real-page-id"


# ---------------------------------------------------------------------------
# _flatten_property — edge cases not covered by the fixture
# ---------------------------------------------------------------------------


def test_flatten_property_empty_dict():
    assert _flatten_property({}) is None


def test_flatten_files_unknown_subtype_silently_skipped():
    prop = {
        "type": "files",
        "files": [
            {"type": "future_cloud", "name": "x", "future_cloud": {"url": "..."}}
        ],
    }
    assert _flatten_property(prop) == []


# ---------------------------------------------------------------------------
# _request — timeout, retry on transient errors, error passthrough
# ---------------------------------------------------------------------------


def test_request_sends_timeout(notion):
    notion._session.request.return_value = _make_response(
        {"results": [], "has_more": False}
    )
    notion.query_database("db")
    assert notion._session.request.call_args[1].get("timeout") == 30


def test_request_retries_502(notion):
    r502 = _make_response({}, status_code=502)
    r200 = _make_response({"results": [], "has_more": False})
    notion._session.request.side_effect = [r502, r200]
    with patch("lamindb.integrations._notion.time") as mock_time:
        mock_time.sleep.return_value = None
        pages = notion.query_database("db")
    assert notion._session.request.call_count == 2
    assert pages == []


def test_request_three_consecutive_429s_exhaust_retries(notion):
    r429 = _make_response({}, status_code=429, headers={"Retry-After": "1"})
    notion._session.request.return_value = r429
    with patch("lamindb.integrations._notion.time") as mock_time:
        mock_time.sleep.return_value = None
        with pytest.raises(RuntimeError, match="rate limit persisted"):
            notion.query_database("db")
    assert notion._session.request.call_count == 3
    r429.raise_for_status.assert_not_called()  # custom message short-circuits


def test_request_retry_after_http_date_falls_back_to_5(notion):
    r429 = _make_response(
        {},
        status_code=429,
        headers={"Retry-After": "Wed, 21 Oct 2024 07:28:00 GMT"},
    )
    r200 = _make_response({"results": [], "has_more": False})
    notion._session.request.side_effect = [r429, r200]
    sleeps: list[int] = []
    with patch("lamindb.integrations._notion.time") as mock_time:
        mock_time.sleep.side_effect = sleeps.append
        notion.query_database("db")
    assert sleeps == [5]


def test_request_500_falls_through_to_raise_for_status(notion):
    r500 = _make_response({}, status_code=500)
    r500.raise_for_status.side_effect = RuntimeError("500 Internal Server Error")
    notion._session.request.return_value = r500
    with pytest.raises(RuntimeError, match="500"):
        notion.query_database("db")
    r500.raise_for_status.assert_called_once()


def test_request_retry_after_large_value_clamped_to_max_backoff(notion):
    r429 = _make_response({}, status_code=429, headers={"Retry-After": "3600"})
    r200 = _make_response({"results": [], "has_more": False})
    notion._session.request.side_effect = [r429, r200]
    sleeps: list[int] = []
    with patch("lamindb.integrations._notion.time") as mock_time:
        mock_time.sleep.side_effect = sleeps.append
        notion.query_database("db")
    assert sleeps == [60]


def test_request_retry_after_negative_clamped_to_zero(notion):
    r429 = _make_response({}, status_code=429, headers={"Retry-After": "-5"})
    r200 = _make_response({"results": [], "has_more": False})
    notion._session.request.side_effect = [r429, r200]
    sleeps: list[int] = []
    with patch("lamindb.integrations._notion.time") as mock_time:
        mock_time.sleep.side_effect = sleeps.append
        notion.query_database("db")
    assert sleeps == [0]
