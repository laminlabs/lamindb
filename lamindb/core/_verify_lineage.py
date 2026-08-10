import ast
import re
from dataclasses import dataclass
from pathlib import Path

# Path matching regex for URIs, absolute/relative paths, and extensions
PATH_PATTERN = re.compile(
    r"^(?:"
    r"[a-zA-Z0-9]+://"                     # URIs (s3://, gs://, https://)
    r"|[a-zA-Z]:[\\/]"                     # Windows drive prefix (C:\, D:/)
    r"|~?[/\\]|\.[\./][/\\]"               # Path prefixes (/, ./, ../, ~/)
    r")[\w\-\.\s/]+$"                      # Path body
    r"|^[\w\-\.\s]+[/\\][\w\-\.\s/]*$"     # Directory relative paths with slashes
    r"|^[\w\-\.\s]+\.[a-zA-Z0-9]{1,10}$"   # Filenames with extensions
)


@dataclass(frozen=True)
class VerifyLineageResult:
    is_fully_tracked: bool
    missing_lineage: tuple[str, ...]


class LaminLineageChecker(ast.NodeVisitor):
    def __init__(self):
        self.tracked_paths: dict[str, list[int]] = {}
        self.untracked_paths: dict[str, list[int]] = {}
        self.var_map: dict[str, set[str]] = {}  # var_name -> set of path strings
        self.artifact_vars: set[str] = set()  # vars assigned from ln.Artifact(...)
        self.has_track_call: bool = False
        self.has_finish_call: bool = False
        self.function_defs: dict[str, ast.FunctionDef] = {}
        self._active_function_calls: set[str] = set()

    def visit_FunctionDef(self, node: ast.FunctionDef):
        self.function_defs[node.name] = node
        self.generic_visit(node)

    def visit_Assign(self, node: ast.Assign):
        """Trace variables assigned to path strings and LaminDB Artifact instances."""
        paths = self._extract_paths_from_node(node.value)
        is_artifact_ctor = (
            isinstance(node.value, ast.Call)
            and self._is_lamindb_artifact_constructor_call(node.value)
        )

        for target in node.targets:
            if isinstance(target, ast.Name):
                if paths:
                    self.var_map[target.id] = paths
                else:
                    self.var_map.pop(target.id, None)
                if is_artifact_ctor:
                    self.artifact_vars.add(target.id)
                else:
                    self.artifact_vars.discard(target.id)
            elif isinstance(target, (ast.Tuple, ast.List)):
                for elt in target.elts:
                    if isinstance(elt, ast.Name):
                        if paths:
                            self.var_map[elt.id] = paths
                        else:
                            self.var_map.pop(elt.id, None)
                        if is_artifact_ctor:
                            self.artifact_vars.add(elt.id)
                        else:
                            self.artifact_vars.discard(elt.id)

        self.generic_visit(node)

    def visit_Call(self, node: ast.Call):
        func_name = self._get_func_name(node.func)
        lineno = getattr(node, "lineno", 0)

        # Record explicit session lifecycle calls used to open/close lineage tracking.
        if func_name in ("ln.track", "lamindb.track"):
            self.has_track_call = True
        if func_name in ("ln.finish", "lamindb.finish"):
            self.has_finish_call = True

        # Distinguish LaminDB API calls from ordinary Python calls.
        # Any path referenced inside LaminDB calls is considered lineage-tracked.
        is_lamin_call = self._is_lamindb_call(node, func_name)

        paths = self._extract_paths_from_node(node)

        if is_lamin_call:
            for path in paths:
                self.tracked_paths.setdefault(path, []).append(lineno)

        else:
            for path in paths:
                self.untracked_paths.setdefault(path, []).append(lineno)

        self._trace_user_function_call(node)
        self.generic_visit(node)

    def _get_func_name(self, node: ast.AST) -> str:
        """Recursively resolves AST attribute names (e.g. ln.Artifact.get)."""
        if isinstance(node, ast.Name):
            return node.id
        elif isinstance(node, ast.Attribute):
            return f"{self._get_func_name(node.value)}.{node.attr}"
        return ""

    def _is_lamindb_artifact_constructor_call(self, node: ast.Call) -> bool:
        func_name = self._get_func_name(node.func)
        return func_name in ("ln.Artifact", "lamindb.Artifact")

    def _is_lamindb_artifact_save_call(self, node: ast.Call) -> bool:
        if not (isinstance(node.func, ast.Attribute) and node.func.attr == "save"):
            return False
        if isinstance(node.func.value, ast.Call):
            return self._is_lamindb_artifact_constructor_call(node.func.value)
        if isinstance(node.func.value, ast.Name):
            return node.func.value.id in self.artifact_vars
        return False

    def _is_lamindb_call(self, node: ast.Call, func_name: str) -> bool:
        if func_name.startswith(("ln.", "lamindb.")):
            return True
        return self._is_lamindb_artifact_save_call(node)

    def _extract_paths_from_node(self, node: ast.AST) -> set[str]:
        """Extracts paths from string constants, variables, and nested function arguments."""
        paths = set()

        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if self._is_path_like(node.value):
                paths.add(node.value)

        elif isinstance(node, ast.Name):
            if node.id in self.var_map:
                paths.update(self.var_map[node.id])

        elif isinstance(node, ast.Call):
            for arg in node.args:
                paths.update(self._extract_paths_from_node(arg))
            for kw in node.keywords:
                paths.update(self._extract_paths_from_node(kw.value))

        return paths

    def _is_path_like(self, s: str) -> bool:
        s = s.strip()
        return bool(s and PATH_PATTERN.match(s))

    def _trace_user_function_call(self, node: ast.Call):
        if not isinstance(node.func, ast.Name):
            return

        function_name = node.func.id
        function_def = self.function_defs.get(function_name)
        if function_def is None or function_name in self._active_function_calls:
            return

        bindings = self._get_param_path_bindings(node, function_def)
        if not bindings:
            return

        self._active_function_calls.add(function_name)
        original_var_map = self.var_map.copy()
        original_artifact_vars = self.artifact_vars.copy()
        try:
            self.var_map.update(bindings)
            for stmt in function_def.body:
                self.visit(stmt)
        finally:
            self.var_map = original_var_map
            self.artifact_vars = original_artifact_vars
            self._active_function_calls.remove(function_name)

    def _get_param_path_bindings(
        self, call: ast.Call, function_def: ast.FunctionDef
    ) -> dict[str, set[str]]:
        bindings: dict[str, set[str]] = {}

        positional_params = function_def.args.posonlyargs + function_def.args.args
        for param, arg in zip(positional_params, call.args):
            paths = self._extract_paths_from_node(arg)
            if paths:
                bindings[param.arg] = paths

        keyword_args = {
            kw.arg: self._extract_paths_from_node(kw.value)
            for kw in call.keywords
            if kw.arg is not None
        }

        for param in positional_params:
            if param.arg in bindings:
                continue
            paths = keyword_args.get(param.arg, set())
            if paths:
                bindings[param.arg] = paths

        for param in function_def.args.kwonlyargs:
            paths = keyword_args.get(param.arg, set())
            if paths:
                bindings[param.arg] = paths

        return bindings


def verify_lineage(script_path: str) -> VerifyLineageResult:
    path = Path(script_path)
    if not path.is_file():
        return VerifyLineageResult(
            is_fully_tracked=False,
            missing_lineage=(f"File not found: {script_path}",),
        )

    with open(path, encoding="utf-8") as f:
        source_code = f.read()

    try:
        tree = ast.parse(source_code, filename=script_path)
    except SyntaxError as e:
        return VerifyLineageResult(
            is_fully_tracked=False,
            missing_lineage=(f"Syntax error while parsing {script_path}: {e}",),
        )

    checker = LaminLineageChecker()
    checker.visit(tree)

    # Keep only paths that are seen in non-LaminDB calls and never seen in LaminDB calls.
    # These are candidates for missing lineage registration.
    truly_untracked_paths = {
        p: lines for p, lines in checker.untracked_paths.items()
        if p not in checker.tracked_paths
    }

    missing_lineage: list[str] = []

    # A script is fully tracked only if it explicitly opens and closes tracking.
    if not checker.has_track_call:
        missing_lineage.append("Missing ln.track() call in script.")

    if not checker.has_finish_call:
        missing_lineage.append("Missing ln.finish() call in script.")

    if truly_untracked_paths:
        for fname, lines in sorted(truly_untracked_paths.items()):
            lines_str = ", ".join(f"line {l}" for l in lines)
            missing_lineage.append(f"File not tracked in lamindb: {fname} ({lines_str})")


    is_fully_tracked = not missing_lineage
    return VerifyLineageResult(
        is_fully_tracked=is_fully_tracked,
        missing_lineage=tuple(missing_lineage),
    )

