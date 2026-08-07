from __future__ import annotations

import ast
import inspect
import os
from dataclasses import dataclass
from pathlib import Path

import lamindb as ln


@dataclass(frozen=True)
class VerifyLineageResult:
    is_fully_tracked: bool
    has_lineage_tracking: bool
    has_external_inputs: bool
    has_external_outputs: bool
    missing_lineage: tuple[str, ...]
    lineage_calls: tuple[str, ...]
    lamindb_input_calls: tuple[str, ...]
    lamindb_output_calls: tuple[str, ...]
    external_input_calls: tuple[str, ...]
    external_output_calls: tuple[str, ...]


def _is_path_like_string(value: str) -> bool:
    if value.startswith(("http://", "https://", "s3://", "gs://")):
        return False
    if "/" in value or "\\" in value:
        return True
    path = Path(value)
    stem, suffix = path.stem, path.suffix
    return bool(stem and suffix and len(suffix) <= 10)


def _normalize_path_token(value: str) -> str:
    return os.path.normpath(value)


def _get_name(node: ast.AST | None) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    return None


def _get_dotted_name(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        left = _get_dotted_name(node.value)
        if left:
            return f"{left}.{node.attr}"
        return node.attr
    return ""


def _extract_call_name(node: ast.Call) -> str:
    return _get_dotted_name(node.func)


def _get_lamindb_output_methods() -> tuple[str, ...]:
    def _get_methods(
        cls: type,
        include: tuple[str, ...] = ("save",),
        prefixes: tuple[str, ...] = ("from_",),
    ) -> list[str]:
        out: list[str] = []
        for name, _obj in inspect.getmembers(cls, predicate=callable):
            if name in include or any(name.startswith(prefix) for prefix in prefixes):
                out.append(name)
        # Some Lamin classes expose callables (notably `save`) that can be
        # callable via getattr but absent from inspect/dir enumeration.
        for name in include:
            if callable(getattr(cls, name, None)):
                out.append(name)
        return sorted(set(out))

    methods: set[str] = set()
    for cls in (ln.Artifact, ln.Collection, ln.Schema):
        methods.update(_get_methods(cls))
    return tuple(sorted(methods))


class _LineageAnalyzer(ast.NodeVisitor):
    def __init__(self) -> None:
        self.lamindb_output_methods = set(_get_lamindb_output_methods())
        self.lamindb_module_aliases: set[str] = set()
        self.imported_symbols: dict[str, str] = {}
        self.var_path_tokens: dict[str, set[str]] = {}
        self.artifact_var_tokens: dict[str, set[str]] = {}

        self.lineage_calls: list[str] = []
        self.lamindb_input_calls: list[str] = []
        self.lamindb_output_calls: list[str] = []
        self.external_input_calls: list[str] = []
        self.external_output_calls: list[str] = []
        self.external_output_tokens: list[set[str]] = []
        self.lamindb_output_tokens: list[set[str]] = []

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            if alias.name == "lamindb":
                self.lamindb_module_aliases.add(alias.asname or alias.name)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        if node.module != "lamindb":
            return
        for alias in node.names:
            self.imported_symbols[alias.asname or alias.name] = alias.name

    def visit_Assign(self, node: ast.Assign) -> None:
        tokens = self._extract_path_tokens(node.value)
        for target in node.targets:
            name = _get_name(target)
            if name is not None:
                if tokens:
                    self.var_path_tokens[name] = set(tokens)
                else:
                    self.var_path_tokens.pop(name, None)
                artifact_tokens = self._extract_artifact_tokens(node.value)
                if artifact_tokens:
                    self.artifact_var_tokens[name] = artifact_tokens
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        call_name = _extract_call_name(node)
        call_text = ast.unparse(node)

        if self._is_lineage_call(node):
            self.lineage_calls.append(call_text)

        if self._is_lamindb_input_call(node):
            self.lamindb_input_calls.append(call_text)

        output_tokens: set[str] = set()
        if self._is_lamindb_output_call(node):
            self.lamindb_output_calls.append(call_text)
            output_tokens = self._extract_output_tokens_from_lamindb_call(node)
            self.lamindb_output_tokens.append(output_tokens)

        if self._is_external_output_call(node):
            self.external_output_calls.append(call_name or call_text)
            write_tokens = self._extract_output_tokens_from_external_call(node)
            self.external_output_tokens.append(write_tokens)

        self.generic_visit(node)

    def _is_lamindb_root(self, node: ast.AST) -> bool:
        if isinstance(node, ast.Name):
            return node.id in self.lamindb_module_aliases
        return False

    def _is_artifact_symbol(self, node: ast.AST) -> bool:
        if isinstance(node, ast.Name):
            return self.imported_symbols.get(node.id) == "Artifact"
        if isinstance(node, ast.Attribute):
            return node.attr == "Artifact" and self._is_lamindb_root(node.value)
        return False

    def _is_track_symbol(self, node: ast.AST) -> bool:
        if isinstance(node, ast.Name):
            return self.imported_symbols.get(node.id) in {"track", "finish"}
        if isinstance(node, ast.Attribute):
            return (
                node.attr in {"track", "finish"}
                and self._is_lamindb_root(node.value)
            )
        return False

    def _is_lineage_call(self, node: ast.Call) -> bool:
        if isinstance(node.func, ast.Name):
            return self.imported_symbols.get(node.func.id) in {"track", "finish"}
        if isinstance(node.func, ast.Attribute):
            return (
                node.func.attr in {"track", "finish"}
                and self._is_lamindb_root(node.func.value)
            )
        return False

    def _is_lamindb_input_call(self, node: ast.Call) -> bool:
        if not isinstance(node.func, ast.Attribute):
            return False
        if node.func.attr not in {"get", "connect"}:
            return False
        return self._is_artifact_symbol(node.func.value)

    def _is_lamindb_output_call(self, node: ast.Call) -> bool:
        if isinstance(node.func, ast.Attribute):
            attr = node.func.attr
            if attr not in self.lamindb_output_methods:
                return False
            if attr.startswith("from_") and self._is_artifact_symbol(node.func.value):
                return True
            if attr == "save":
                if isinstance(node.func.value, ast.Call):
                    return self._is_artifact_symbol(node.func.value.func)
                if isinstance(node.func.value, ast.Name):
                    return node.func.value.id in self.artifact_var_tokens
        return False

    def _is_external_output_call(self, node: ast.Call) -> bool:
        if isinstance(node.func, ast.Name) and node.func.id == "open":
            if len(node.args) >= 2 and self._is_write_mode(node.args[1]):
                return True
            for kw in node.keywords:
                if kw.arg == "mode" and self._is_write_mode(kw.value):
                    return True

        if isinstance(node.func, ast.Attribute):
            if node.func.attr in {"write_text", "write_bytes", "touch"}:
                return True
            dotted = _get_dotted_name(node.func)
            if dotted.endswith(".save"):
                if isinstance(node.func.value, ast.Name):
                    if node.func.value.id in self.artifact_var_tokens:
                        return False
                if isinstance(node.func.value, ast.Call):
                    if self._is_artifact_symbol(node.func.value.func):
                        return False
                return True
            if dotted.endswith(".savetxt"):
                return True
            if dotted.endswith(".savez") or dotted.endswith(".savez_compressed"):
                return True
        return False

    def _is_write_mode(self, node: ast.AST) -> bool:
        if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
            return False
        mode = node.value
        return any(flag in mode for flag in ("w", "a", "x", "+"))

    def _extract_artifact_tokens(self, node: ast.AST) -> set[str]:
        if not isinstance(node, ast.Call):
            return set()
        if not self._is_artifact_symbol(node.func):
            return set()
        tokens: set[str] = set()
        if node.args:
            tokens.update(self._extract_path_tokens(node.args[0]))
        for kw in node.keywords:
            if kw.arg == "key":
                tokens.update(self._extract_path_tokens(kw.value))
        return tokens

    def _extract_output_tokens_from_lamindb_call(self, node: ast.Call) -> set[str]:
        if isinstance(node.func, ast.Attribute) and node.func.attr == "save":
            if isinstance(node.func.value, ast.Call):
                return self._extract_artifact_tokens(node.func.value)
            if isinstance(node.func.value, ast.Name):
                return set(self.artifact_var_tokens.get(node.func.value.id, set()))

        if (
            isinstance(node.func, ast.Attribute)
            and node.func.attr.startswith("from_")
            and self._is_artifact_symbol(node.func.value)
        ):
            tokens: set[str] = set()
            for kw in node.keywords:
                if kw.arg == "key":
                    tokens.update(self._extract_path_tokens(kw.value))
            return tokens
        return set()

    def _extract_output_tokens_from_external_call(self, node: ast.Call) -> set[str]:
        if isinstance(node.func, ast.Name) and node.func.id == "open" and node.args:
            return self._extract_path_tokens(node.args[0])
        if isinstance(node.func, ast.Attribute):
            attr = node.func.attr
            if attr in {"write_text", "write_bytes", "touch"}:
                return self._extract_path_tokens(node.func.value)
            if attr in {"save", "savetxt", "savez", "savez_compressed"} and node.args:
                return self._extract_path_tokens(node.args[0])
        return set()

    def _extract_path_tokens(self, node: ast.AST) -> set[str]:
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if _is_path_like_string(node.value):
                return {_normalize_path_token(node.value)}
            return set()
        if isinstance(node, ast.Name):
            if node.id in self.var_path_tokens:
                return set(self.var_path_tokens[node.id])
            return {f"$var:{node.id}"}
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id in {
                "Path",
                "PurePath",
                "PosixPath",
                "WindowsPath",
            }:
                tokens: set[str] = set()
                for arg in node.args:
                    tokens.update(self._extract_path_tokens(arg))
                return tokens
            if isinstance(node.func, ast.Attribute) and node.func.attr == "join":
                tokens: set[str] = set()
                for arg in node.args:
                    tokens.update(self._extract_path_tokens(arg))
                return tokens
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return self._extract_path_tokens(node.left) | self._extract_path_tokens(
                node.right
            )
        return set()


def verify_lineage(script_path: str | Path) -> VerifyLineageResult:
    path = Path(script_path)
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))

    analyzer = _LineageAnalyzer()
    analyzer.visit(tree)

    has_track = any("track" in call for call in analyzer.lineage_calls)
    has_finish = any("finish" in call for call in analyzer.lineage_calls)
    has_lineage_tracking = has_track and has_finish

    lamindb_output_union: set[str] = set()
    for tokens in analyzer.lamindb_output_tokens:
        lamindb_output_union.update(tokens)

    untracked_output_calls: list[str] = []
    for idx, call in enumerate(analyzer.external_output_calls):
        write_tokens = analyzer.external_output_tokens[idx]
        if write_tokens and write_tokens & lamindb_output_union:
            continue
        untracked_output_calls.append(call)

    has_external_outputs = len(untracked_output_calls) > 0
    external_output_calls = tuple(untracked_output_calls)

    missing_lineage: list[str] = []
    if not has_lineage_tracking:
        missing_lineage.append("missing lineage tracking call: expected track() and finish()")
    if has_external_outputs:
        missing_lineage.append("unexpected non-LaminDB output writes detected")

    is_fully_tracked = has_lineage_tracking and not has_external_outputs

    return VerifyLineageResult(
        is_fully_tracked=is_fully_tracked,
        has_lineage_tracking=has_lineage_tracking,
        has_external_inputs=False,
        has_external_outputs=has_external_outputs,
        missing_lineage=tuple(missing_lineage),
        lineage_calls=tuple(analyzer.lineage_calls),
        lamindb_input_calls=tuple(analyzer.lamindb_input_calls),
        lamindb_output_calls=tuple(analyzer.lamindb_output_calls),
        external_input_calls=tuple(analyzer.external_input_calls),
        external_output_calls=external_output_calls,
    )
