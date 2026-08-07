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

EXTERNAL_OUTPUT_PATH_ARG_POSITIONS = {
    "numpy.save": 0,
    "np.save": 0,
    "numpy.savetxt": 0,
    "np.savetxt": 0,
    "matplotlib.pyplot.savefig": 0,
    "plt.savefig": 0,
}

EXTERNAL_OUTPUT_PATH_ARG_KEYWORDS = {
    "numpy.save": {"file"},
    "np.save": {"file"},
    "numpy.savetxt": {"fname"},
    "np.savetxt": {"fname"},
    "matplotlib.pyplot.savefig": {"fname"},
    "plt.savefig": {"fname"},
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
    missing_lineage: tuple[str, ...]
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


def _saved_artifact_arg_lines(
    tree: ast.AST, lamindb_module_aliases: set[str], lamindb_model_aliases: set[str]
) -> dict[str, list[int]]:
    saved_arg_lines: dict[str, list[int]] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute) or node.func.attr not in LAMINDB_OUTPUT_METHODS:
            continue

        artifact_constructor_call = node.func.value
        if not isinstance(artifact_constructor_call, ast.Call):
            continue
        constructor_call_name = _dotted_name(artifact_constructor_call.func)
        if constructor_call_name is None:
            continue

        root = _root_name(artifact_constructor_call.func)
        is_lamindb_artifact_constructor = False
        if root in lamindb_model_aliases and constructor_call_name.split(".")[-1] == "Artifact":
            is_lamindb_artifact_constructor = True
        elif root in lamindb_module_aliases:
            call_parts = constructor_call_name.split(".")
            if len(call_parts) > 1 and call_parts[1] == "Artifact":
                is_lamindb_artifact_constructor = True
        if not is_lamindb_artifact_constructor or len(artifact_constructor_call.args) == 0:
            continue

        first_arg = artifact_constructor_call.args[0]
        if isinstance(first_arg, ast.Name):
            saved_arg_lines.setdefault(first_arg.id, []).append(node.lineno)
    return saved_arg_lines


def _is_later_saved_path_write(node: ast.Call, saved_artifact_arg_lines: dict[str, list[int]]) -> bool:
    if not isinstance(node.func, ast.Attribute):
        return False
    if node.func.attr not in {"write_text", "write_bytes"}:
        return False
    if not isinstance(node.func.value, ast.Name):
        return False

    candidate_saved_lines = saved_artifact_arg_lines.get(node.func.value.id, [])
    return any(save_line > node.lineno for save_line in candidate_saved_lines)


def _saved_name_argument(
    node: ast.Call, *, arg_position: int, keyword_names: set[str]
) -> str | None:
    if len(node.args) > arg_position and isinstance(node.args[arg_position], ast.Name):
        return node.args[arg_position].id
    for keyword in node.keywords:
        if keyword.arg in keyword_names and isinstance(keyword.value, ast.Name):
            return keyword.value.id
    return None


def _is_later_saved_external_output_call(
    node: ast.Call, call_name: str, saved_artifact_arg_lines: dict[str, list[int]]
) -> bool:
    arg_position = EXTERNAL_OUTPUT_PATH_ARG_POSITIONS.get(call_name)
    if arg_position is None:
        return False
    path_name = _saved_name_argument(
        node,
        arg_position=arg_position,
        keyword_names=EXTERNAL_OUTPUT_PATH_ARG_KEYWORDS.get(call_name, set()),
    )
    if path_name is None:
        return False
    candidate_saved_lines = saved_artifact_arg_lines.get(path_name, [])
    return any(save_line > node.lineno for save_line in candidate_saved_lines)


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

    saved_artifact_arg_lines = _saved_artifact_arg_lines(
        tree=tree,
        lamindb_module_aliases=lamindb_module_aliases,
        lamindb_model_aliases=lamindb_model_aliases,
    )

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
            if _is_later_saved_external_output_call(
                node, call_name, saved_artifact_arg_lines
            ):
                continue
            external_output_calls.append(_format_call(call_name, node.lineno))
            continue
        if method_name in EXTERNAL_OUTPUT_METHOD_NAMES:
            if _is_later_saved_path_write(node, saved_artifact_arg_lines):
                continue
            external_output_calls.append(_format_call(call_name, node.lineno))

    has_lineage_tracking = len(lineage_calls) > 0
    has_external_inputs = len(external_input_calls) > 0
    has_external_outputs = len(external_output_calls) > 0

    missing_lineage: list[str] = []
    if not has_lineage_tracking:
        missing_lineage.append("lineage tracking call (`ln.track()` or `ln.finish()`)")
    if has_external_inputs:
        missing_lineage.append(
            _summarize_calls(
                "unexpected non-LaminDB input reads detected",
                external_input_calls,
            )
        )
    if has_external_outputs:
        missing_lineage.append(
            _summarize_calls(
                "unexpected non-LaminDB output writes detected",
                external_output_calls,
            )
        )

    return ScriptLineageVerification(
        has_lineage_tracking=has_lineage_tracking,
        has_external_inputs=has_external_inputs,
        has_external_outputs=has_external_outputs,
        is_fully_tracked=len(missing_lineage) == 0,
        missing_lineage=tuple(missing_lineage),
        lineage_calls=tuple(lineage_calls),
        lamindb_input_calls=tuple(lamindb_input_calls),
        lamindb_output_calls=tuple(lamindb_output_calls),
        external_input_calls=tuple(external_input_calls),
        external_output_calls=tuple(external_output_calls),
    )
