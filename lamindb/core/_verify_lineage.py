from __future__ import annotations

import ast
import pathlib
from dataclasses import dataclass


LAMINDB_MODEL_NAMES = {
    "Artifact",
    "Collection",
    "Run",
    "Transform",
}

LAMINDB_INPUT_METHODS = {
    "get",
    "filter",
    "search",
    "one",
    "one_or_none",
    "first",
    "last",
    "load",
    "open",
    "cache",
    "mapped",
}

LAMINDB_OUTPUT_METHODS = {
    "save",
}

EXTERNAL_INPUT_CALL_NAMES = {
    "open",
    "numpy.load",
    "np.load",
    "anndata.read",
    "ad.read",
    "scanpy.read",
    "sc.read",
    "sqlite3.connect",
    "duckdb.connect",
    "fsspec.open",
}

EXTERNAL_INPUT_PREFIXES = (
    "pandas.read_",
    "pd.read_",
    "polars.read_",
    "pl.read_",
)

EXTERNAL_OUTPUT_CALL_NAMES = {
    "numpy.save",
    "np.save",
    "numpy.savetxt",
    "np.savetxt",
    "pickle.dump",
    "json.dump",
    "yaml.dump",
    "toml.dump",
    "matplotlib.pyplot.savefig",
    "plt.savefig",
}

EXTERNAL_OUTPUT_METHOD_NAMES = {
    "to_csv",
    "to_parquet",
    "to_json",
    "to_excel",
    "to_pickle",
    "to_feather",
    "to_sql",
    "to_hdf",
    "savefig",
    "write_text",
    "write_bytes",
    "write",
    "writelines",
}


@dataclass(frozen=True)
class ScriptLineageVerification:
    """Result of static lineage checks for a Python script."""

    has_lineage_tracking: bool
    has_external_inputs: bool
    has_external_outputs: bool
    is_fully_tracked: bool
    missing: tuple[str, ...]
    lineage_calls: tuple[str, ...]
    lamindb_input_calls: tuple[str, ...]
    lamindb_output_calls: tuple[str, ...]
    external_input_calls: tuple[str, ...]
    external_output_calls: tuple[str, ...]


def _dotted_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _dotted_name(node.value)
        if parent is None:
            return None
        return f"{parent}.{node.attr}"
    if isinstance(node, ast.Call):
        return _dotted_name(node.func)
    return None


def _root_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return _root_name(node.value)
    if isinstance(node, ast.Call):
        return _root_name(node.func)
    return None


def _format_call(call_name: str, lineno: int) -> str:
    return f"{call_name} (line {lineno})"


def _is_open_output_mode(node: ast.Call) -> bool:
    mode_value: str | None = None
    if len(node.args) > 1 and isinstance(node.args[1], ast.Constant):
        if isinstance(node.args[1].value, str):
            mode_value = node.args[1].value
    for keyword in node.keywords:
        if keyword.arg == "mode" and isinstance(keyword.value, ast.Constant):
            if isinstance(keyword.value.value, str):
                mode_value = keyword.value.value
                break
    if mode_value is None:
        mode_value = "r"
    return any(flag in mode_value for flag in ("w", "a", "x", "+"))


def _summarize_calls(prefix: str, calls: list[str]) -> str:
    return f"{prefix}: {', '.join(calls)}"


def verify_lineage(path: str | pathlib.Path) -> ScriptLineageVerification:
    """Statically verify lineage tracking conventions in a Python script.

    Checks for:
    - run lineage calls (`ln.track(...)` and `ln.finish(...)`)
    - laminDB input retrieval calls (for example `ln.Artifact.get(...)`)
    - laminDB output persistence calls (for example `.save()` on laminDB records)
    - common non-laminDB input reads (for example `pd.read_csv(...)`)
    - common non-laminDB output writes (for example `open(..., "w")`)

    Notes:
    - This is a static AST check and cannot prove runtime behavior.
    - Zero-input and zero-output scripts are valid and can pass.
    - "All I/O comes from LaminDB" is interpreted conservatively:
      if known non-LaminDB read/write patterns are present, the check fails.
    """
    script_path = pathlib.Path(path)
    source = script_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    lamindb_module_aliases: set[str] = set()
    lamindb_model_aliases: set[str] = set()
    imported_track_aliases: set[str] = set()
    imported_finish_aliases: set[str] = set()
    imported_save_aliases: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "lamindb":
                    lamindb_module_aliases.add(alias.asname or alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module == "lamindb":
                for alias in node.names:
                    imported_name = alias.asname or alias.name
                    if alias.name == "track":
                        imported_track_aliases.add(imported_name)
                    elif alias.name == "finish":
                        imported_finish_aliases.add(imported_name)
                    elif alias.name == "save":
                        imported_save_aliases.add(imported_name)
                    elif alias.name in LAMINDB_MODEL_NAMES:
                        lamindb_model_aliases.add(imported_name)
            elif node.module == "lamindb.models":
                for alias in node.names:
                    if alias.name in LAMINDB_MODEL_NAMES:
                        lamindb_model_aliases.add(alias.asname or alias.name)

    lineage_calls: list[str] = []
    lamindb_input_calls: list[str] = []
    lamindb_output_calls: list[str] = []
    external_input_calls: list[str] = []
    external_output_calls: list[str] = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue

        call_name = _dotted_name(node.func)
        if call_name is None:
            continue

        # Lineage tracking calls
        if (
            call_name in imported_track_aliases
            or any(call_name == f"{alias}.track" for alias in lamindb_module_aliases)
            or call_name in imported_finish_aliases
            or any(call_name == f"{alias}.finish" for alias in lamindb_module_aliases)
        ):
            lineage_calls.append(_format_call(call_name, node.lineno))

        # LaminDB inputs
        method_name = call_name.split(".")[-1]
        if method_name in LAMINDB_INPUT_METHODS:
            root = _root_name(node.func)
            if root in lamindb_model_aliases:
                lamindb_input_calls.append(_format_call(call_name, node.lineno))
            elif root in lamindb_module_aliases:
                # Handles `ln.Artifact.get(...)`, `ln.Collection.filter(...)`, etc.
                if "." in call_name:
                    model_name = call_name.split(".")[1]
                    if model_name in LAMINDB_MODEL_NAMES:
                        lamindb_input_calls.append(_format_call(call_name, node.lineno))

        # LaminDB outputs
        if call_name in imported_save_aliases:
            lamindb_output_calls.append(_format_call(call_name, node.lineno))
        elif method_name in LAMINDB_OUTPUT_METHODS:
            root = _root_name(node.func)
            if root in lamindb_model_aliases:
                lamindb_output_calls.append(_format_call(call_name, node.lineno))
            elif root in lamindb_module_aliases:
                # Handles `ln.Artifact(...).save()` and related calls.
                lamindb_output_calls.append(_format_call(call_name, node.lineno))

        # External inputs and outputs
        if call_name in {"open", "fsspec.open"}:
            formatted = _format_call(call_name, node.lineno)
            if _is_open_output_mode(node):
                external_output_calls.append(formatted)
            else:
                external_input_calls.append(formatted)
            continue

        if call_name in EXTERNAL_INPUT_CALL_NAMES:
            external_input_calls.append(_format_call(call_name, node.lineno))
            continue
        if call_name.endswith(".read_text") or call_name.endswith(".read_bytes"):
            external_input_calls.append(_format_call(call_name, node.lineno))
            continue
        if call_name.startswith(EXTERNAL_INPUT_PREFIXES):
            external_input_calls.append(_format_call(call_name, node.lineno))
            continue

        if call_name in EXTERNAL_OUTPUT_CALL_NAMES:
            external_output_calls.append(_format_call(call_name, node.lineno))
            continue
        if method_name in EXTERNAL_OUTPUT_METHOD_NAMES:
            external_output_calls.append(_format_call(call_name, node.lineno))

    has_lineage_tracking = len(lineage_calls) > 0
    has_external_inputs = len(external_input_calls) > 0
    has_external_outputs = len(external_output_calls) > 0

    missing: list[str] = []
    if not has_lineage_tracking:
        missing.append("lineage tracking call (`ln.track()` or `ln.finish()`)")
    if has_external_inputs:
        missing.append(
            _summarize_calls(
                "unexpected non-LaminDB input reads detected",
                external_input_calls,
            )
        )
    if has_external_outputs:
        missing.append(
            _summarize_calls(
                "unexpected non-LaminDB output writes detected",
                external_output_calls,
            )
        )

    return ScriptLineageVerification(
        has_lineage_tracking=has_lineage_tracking,
        has_external_inputs=has_external_inputs,
        has_external_outputs=has_external_outputs,
        is_fully_tracked=len(missing) == 0,
        missing=tuple(missing),
        lineage_calls=tuple(lineage_calls),
        lamindb_input_calls=tuple(lamindb_input_calls),
        lamindb_output_calls=tuple(lamindb_output_calls),
        external_input_calls=tuple(external_input_calls),
        external_output_calls=tuple(external_output_calls),
    )
