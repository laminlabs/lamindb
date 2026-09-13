from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Literal

import lamindb as ln
from lamindb.core._context import (
    REDACTED_SECRET_VALUE,
    _validate_with_virtual_schema,
    serialize_params_to_json,
)
from lamindb_setup.core.upath import UPath


def test_serialize_params_to_json():
    a_path = Path("/some/local/folder")
    a_upath = UPath("s3://bucket/key")
    params = {
        "path_key": a_path,
        "none_key": None,
        "empty_list_key": [],
        "list_str_key": ["string"],
        "upath_key": a_upath,
        "str_key": "plain",
        "api_key": "test-api-key-value",
        "openAIApiKey": "another-secret",
        "database_url": "postgresql://db_user:db_password@db.example.com:5432/mydb",
    }
    result = serialize_params_to_json(params)
    # None is omitted
    assert "none_key" not in result
    # Empty list is omitted (same as None)
    assert "empty_list_key" not in result
    # Path is serialized to posix string
    assert result["path_key"] == "/some/local/folder"
    # UPath is serialized to posix string
    assert result["upath_key"] == "s3://bucket/key"
    # List of strings is JSON-serialized as-is (list[cat ? str])
    assert result["list_str_key"] == ["string"]
    # Other values unchanged
    assert result["str_key"] == "plain"
    assert result["api_key"] == REDACTED_SECRET_VALUE
    assert result["openAIApiKey"] == REDACTED_SECRET_VALUE
    assert result["database_url"] == REDACTED_SECRET_VALUE
    assert set(result.keys()) == {
        "path_key",
        "upath_key",
        "str_key",
        "list_str_key",
        "api_key",
        "openAIApiKey",
        "database_url",
    }


def test_serialize_params_to_json_redacts_provider_api_key_names():
    params = {
        "LAMIN_API_KEY": "lamin-super-secret",
        "OPENAI_API_KEY": "openai-super-secret",
        "ANTHROPIC_API_KEY": "anthropic-super-secret",
        "GEMINI_API_KEY": "gemini-super-secret",
        "provider_name": "safe-value",
    }
    result = serialize_params_to_json(params)
    assert result["LAMIN_API_KEY"] == REDACTED_SECRET_VALUE
    assert result["OPENAI_API_KEY"] == REDACTED_SECRET_VALUE
    assert result["ANTHROPIC_API_KEY"] == REDACTED_SECRET_VALUE
    assert result["GEMINI_API_KEY"] == REDACTED_SECRET_VALUE
    assert result["provider_name"] == "safe-value"


def test_serialize_params_to_json_skips_annotation_mismatch(ccaplog):
    params = {"count": "not-an-int", "label": "ok"}
    result = serialize_params_to_json(
        params, expected_param_types={"count": int, "label": str}
    )
    assert result == {"label": "ok"}
    assert "does not match annotation" in ccaplog.text
    assert "count" in ccaplog.text


def test_serialize_params_to_json_serializes_valid_record_annotations():
    record = ln.Record(name="track-param-record").save()
    try:
        params = {"record": record, "records": [record]}
        result = serialize_params_to_json(
            params,
            expected_param_types={"record": ln.Record, "records": list[ln.Record]},
        )
        assert result == {
            "record": f"Record[{record.uid}]",
            "records": [f"Record[{record.uid}]"],
        }
    finally:
        record.delete(permanent=True)


def test_serialize_params_to_json_skips_unsupported_annotation(ccaplog):
    params = {"mode": "turbo"}
    result = serialize_params_to_json(
        params, expected_param_types={"mode": Literal["fast", "slow"]}
    )
    assert "mode" not in result
    assert "literal value mismatch" in ccaplog.text


def test_serialize_params_to_json_supports_literal_annotation():
    params = {"mode": "fast"}
    result = serialize_params_to_json(
        params, expected_param_types={"mode": Literal["fast", "slow"]}
    )
    assert result == {"mode": "fast"}


def test_serialize_params_to_json_supports_sequence_mapping_and_nested_union():
    params = {
        "names": ["a", "b"],
        "counts": {"a": 1, "b": 2},
        "mixed": [1, "two", 3],
    }
    result = serialize_params_to_json(
        params,
        expected_param_types={
            "names": Sequence[str],
            "counts": Mapping[str, int],
            "mixed": list[int | str],
        },
    )
    assert result == params


def test_serialize_params_to_json_skips_unresolved_string_annotation(ccaplog):
    params = {"value": "alpha"}
    result = serialize_params_to_json(
        params, expected_param_types={"value": "UnknownType"}
    )
    assert result == {}
    assert "unresolved string annotation" in ccaplog.text


def test_virtual_schema_validation_uses_params_key_label(monkeypatch):
    captured: dict[str, str] = {}

    class DummyCurator:
        def __init__(
            self,
            dataset,
            schema,
            slot=None,
            require_saved_schema=False,
            using=None,
            key_label="items",
        ):
            captured["key_label"] = key_label

        def validate(self):
            return None

    monkeypatch.setattr("lamindb.curators.core.ExperimentalDictCurator", DummyCurator)
    valid_keys, invalid_reasons = _validate_with_virtual_schema(
        {"count": 1}, {"count": int}
    )

    assert captured["key_label"] == "params"
    assert valid_keys == {"count"}
    assert invalid_reasons == {}
