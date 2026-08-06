from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path


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


@dataclass(frozen=True)
class ScriptLineageVerification:
    """Result of static lineage checks for a Python script."""

    has_lineage_tracking: bool
    has_lamindb_inputs: bool
    has_lamindb_outputs: bool
    has_external_inputs: bool
    is_fully_tracked: bool
    missing: tuple[str, ...]
    lineage_calls: tuple[str, ...]
    lamindb_input_calls: tuple[str, ...]
    lamindb_output_calls: tuple[str, ...]
    external_input_calls: tuple[str, ...]


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


def verify_lineage(path: str | Path) -> ScriptLineageVerification:
    """Statically verify lineage tracking conventions in a Python script.

    Checks for:
    - run lineage calls (`ln.track(...)` and `ln.finish(...)`)
    - laminDB input retrieval calls (for example `ln.Artifact.get(...)`)
    - laminDB output persistence calls (for example `.save()` on laminDB records)
    - common non-laminDB input reads (for example `pd.read_csv(...)`)

    Notes:
    - This is a static AST check and cannot prove runtime behavior.
    - "All inputs come from LaminDB" is interpreted conservatively:
      if known external read patterns are present, the check fails.
    """
    script_path = Path(path)
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

        # External inputs
        if call_name in EXTERNAL_INPUT_CALL_NAMES:
            external_input_calls.append(_format_call(call_name, node.lineno))
            continue
        if call_name.endswith(".read_text") or call_name.endswith(".read_bytes"):
            external_input_calls.append(_format_call(call_name, node.lineno))
            continue
        if call_name.startswith(EXTERNAL_INPUT_PREFIXES):
            external_input_calls.append(_format_call(call_name, node.lineno))

    has_lineage_tracking = len(lineage_calls) > 0
    has_lamindb_inputs = len(lamindb_input_calls) > 0
    has_lamindb_outputs = len(lamindb_output_calls) > 0
    has_external_inputs = len(external_input_calls) > 0

    missing: list[str] = []
    if not has_lineage_tracking:
        missing.append("lineage tracking call (`ln.track()` or `ln.finish()`)")
    if not has_lamindb_inputs:
        missing.append("LaminDB input retrieval")
    if not has_lamindb_outputs:
        missing.append("LaminDB output persistence")
    if has_external_inputs:
        missing.append("unexpected non-LaminDB input reads")

    return ScriptLineageVerification(
        has_lineage_tracking=has_lineage_tracking,
        has_lamindb_inputs=has_lamindb_inputs,
        has_lamindb_outputs=has_lamindb_outputs,
        has_external_inputs=has_external_inputs,
        is_fully_tracked=len(missing) == 0,
        missing=tuple(missing),
        lineage_calls=tuple(lineage_calls),
        lamindb_input_calls=tuple(lamindb_input_calls),
        lamindb_output_calls=tuple(lamindb_output_calls),
        external_input_calls=tuple(external_input_calls),
    )
