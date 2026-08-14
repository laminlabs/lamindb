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

# Calls in these categories can carry path strings but are considered
# environment/path setup (not data lineage I/O), so they are ignored.
NON_LINEAGE_EXACT_CALLS = frozenset(
    {
        "Path",
        "PurePath",
        "PosixPath",
        "WindowsPath",
        "makedirs",
        "os.makedirs",
    }
)
NON_LINEAGE_CALL_PREFIXES = (
    "sys.path.",
    "os.path.",
    "importlib.",
)
NON_LINEAGE_ATTR_CALLS = frozenset(
    {
        "Path",
        "PurePath",
        "PosixPath",
        "WindowsPath",
        "mkdir",
        "exists",
        "is_file",
        "is_dir",
        "resolve",
        "expanduser",
        "glob",
        "rglob",
    }
)


@dataclass(frozen=True)
class VerifyLineageResult:
    """Public result object returned by `verify_lineage()`."""
    is_fully_tracked: bool
    missing_lineage: tuple[str, ...]


class LaminLineageChecker(ast.NodeVisitor):
    def __init__(self):
        # Paths observed in LaminDB API calls, keyed by path string -> line numbers.
        self.tracked_paths: dict[str, list[int]] = {}
        # Paths observed in non-LaminDB calls, keyed by path string -> line numbers.
        self.untracked_paths: dict[str, list[int]] = {}
        # Best-effort binding of variable names to path-like strings as code is traversed.
        self.var_map: dict[str, set[str]] = {}  # var_name -> set of path strings
        # Variable names that currently refer to an ln.Artifact(...) instance.
        self.artifact_vars: set[str] = set()  # vars assigned from ln.Artifact(...)
        # Session-level lifecycle flags required for a "fully tracked" script.
        self.has_track_call: bool = False
        self.has_finish_call: bool = False
        # Module aliases bound to `lamindb` via `import lamindb as <alias>`.
        self.lamindb_module_aliases: set[str] = {"lamindb"}
        # Local symbols imported from `lamindb`, keyed by local alias -> imported name.
        self.imported_lamindb_symbols: dict[str, str] = {}
        # Function definitions are cached so later calls can be re-visited with bound args.
        self.function_defs: dict[str, ast.FunctionDef] = {}
        # Prevent recursive re-entry while tracing user-defined function calls.
        self._active_function_calls: set[str] = set()

    def visit_ImportFrom(self, node: ast.ImportFrom):
        if node.module != "lamindb":
            return

        for imported in node.names:
            local_name = imported.asname or imported.name
            self.imported_lamindb_symbols[local_name] = imported.name

    def visit_Import(self, node: ast.Import):
        for imported in node.names:
            if imported.name == "lamindb":
                local_name = imported.asname or imported.name
                self.lamindb_module_aliases.add(local_name)

    def visit_FunctionDef(self, node: ast.FunctionDef):
        for decorator in node.decorator_list:
            decorator_name = (
                self._get_func_name(decorator.func)
                if isinstance(decorator, ast.Call)
                else self._get_func_name(decorator)
            )
            if self._is_lamindb_flow_name(decorator_name):
                # `@ln.flow` provides run lifecycle management.
                self.has_track_call = True
                self.has_finish_call = True
        # Save function bodies for later path-flow tracing at call sites.
        self.function_defs[node.name] = node
        self.generic_visit(node)

    def visit_Assign(self, node: ast.Assign):
        """Trace variables assigned to path strings and LaminDB Artifact instances."""
        if self._is_env_var_setup_assignment(node):
            # Environment variable assignment is runtime setup and not lineage I/O.
            # Skip traversing nested helper calls (e.g. str(path)) in this subtree.
            return

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

    def _is_env_var_setup_assignment(self, node: ast.Assign) -> bool:
        for target in node.targets:
            if not isinstance(target, ast.Subscript):
                continue
            if self._get_func_name(target.value) == "os.environ":
                return True
        return False

    def visit_Call(self, node: ast.Call):
        func_name = self._get_func_name(node.func)
        lineno = getattr(node, "lineno", 0)

        # Record explicit session lifecycle calls used to open/close lineage tracking.
        if self._is_lamindb_module_call(func_name, "track") or self._is_imported_lamindb_symbol(
            func_name, "track"
        ):
            self.has_track_call = True
        if self._is_lamindb_module_call(func_name, "finish") or self._is_imported_lamindb_symbol(
            func_name, "finish"
        ):
            self.has_finish_call = True
        if self._is_lamindb_flow_name(func_name):
            self.has_track_call = True
            self.has_finish_call = True

        # Distinguish LaminDB API calls from ordinary Python calls.
        # Any path referenced inside LaminDB calls is considered lineage-tracked.
        is_lamin_call = self._is_lamindb_call(node, func_name)
        is_non_lineage_setup_call = self._is_non_lineage_setup_call(node, func_name)
        if is_non_lineage_setup_call:
            # Ignore full setup-call subtrees (e.g. `sys.path.insert(..., str(path))`)
            # so nested helper calls are not misclassified as lineage-relevant I/O.
            return

        paths = self._extract_paths_from_node(node)

        if is_lamin_call:
            for path in paths:
                self.tracked_paths.setdefault(path, []).append(lineno)

        elif not is_non_lineage_setup_call:
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
        return self._is_lamindb_module_call(func_name, "Artifact") or self._is_imported_lamindb_symbol(
            func_name, "Artifact"
        )

    def _is_lamindb_module_call(self, func_name: str, call_name: str) -> bool:
        root_name, _, remainder = func_name.partition(".")
        return root_name in self.lamindb_module_aliases and remainder == call_name

    def _is_imported_lamindb_symbol(self, func_name: str, symbol_name: str) -> bool:
        root_name = func_name.split(".", 1)[0]
        return self.imported_lamindb_symbols.get(root_name) == symbol_name

    def _is_lamindb_flow_name(self, func_name: str) -> bool:
        return self._is_lamindb_module_call(func_name, "flow") or self._is_imported_lamindb_symbol(
            func_name, "flow"
        )

    def _is_lamindb_artifact_save_call(self, node: ast.Call) -> bool:
        # Treat both `ln.Artifact(...).save()` and `artifact.save()` (where artifact
        # was created from ln.Artifact(...)) as LaminDB lineage-tracked operations.
        if not (isinstance(node.func, ast.Attribute) and node.func.attr == "save"):
            return False
        if isinstance(node.func.value, ast.Call):
            return self._is_lamindb_artifact_constructor_call(node.func.value)
        if isinstance(node.func.value, ast.Name):
            return node.func.value.id in self.artifact_vars
        return False

    def _is_lamindb_call(self, node: ast.Call, func_name: str) -> bool:
        # All ln./lamindb.* calls are considered lineage-aware.
        # We also include artifact save calls recognized above.
        root_name = func_name.split(".", 1)[0]
        if root_name in self.lamindb_module_aliases:
            return True
        if root_name in self.imported_lamindb_symbols:
            return True
        return self._is_lamindb_artifact_save_call(node)

    def _is_non_lineage_setup_call(self, node: ast.Call, func_name: str) -> bool:
        # Centralized allowlist for path/setup helpers that should not be
        # interpreted as lineage-relevant file I/O.
        if func_name in NON_LINEAGE_EXACT_CALLS:
            return True
        if func_name.startswith(NON_LINEAGE_CALL_PREFIXES):
            return True
        return isinstance(node.func, ast.Attribute) and node.func.attr in NON_LINEAGE_ATTR_CALLS

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

        elif isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Div)):
            # Handle both string concatenation and pathlib-style joins.
            left_paths = self._extract_paths_from_node(node.left)
            right_paths = self._extract_paths_from_node(node.right)
            for left_path in left_paths:
                for right_path in right_paths:
                    combined = self._join_paths(left_path, right_path)
                    if self._is_path_like(combined):
                        paths.add(combined)

        return paths

    def _is_path_like(self, s: str) -> bool:
        s = s.strip()
        return bool(s and PATH_PATTERN.match(s))

    def _join_paths(self, left: str, right: str) -> str:
        left = left.strip()
        right = right.strip()
        if left.endswith(("/", "\\")) or right.startswith(("/", "\\")):
            return f"{left}{right}"
        return f"{left}/{right}"

    def _trace_user_function_call(self, node: ast.Call):
        """Propagate path bindings into user-defined functions at call sites.

        This gives the checker a lightweight inter-procedural view, so path usage
        hidden behind helper functions can still be classified as tracked/untracked.
        """
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
        # Temporarily augment variable scope with call argument bindings, then
        # restore state after visiting the function body.
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
    """Statically analyze one script and report missing LaminDB lineage hooks.

    Coverage (static AST checks):
    - Confirms the script contains explicit session lifecycle calls:
      `ln.track()` and `ln.finish()`.
    - Extracts path-like strings from literals, variables, call arguments,
      and simple path joins (string `+` or path `/` operations).
    - Classifies path usage as tracked when it appears in LaminDB calls
      (`ln.*`, `lamindb.*`, including `ln.Artifact(...).save()` patterns),
      and flags paths that only appear in non-LaminDB calls.
    - Propagates path-like argument bindings into user-defined helper
      functions to catch indirect path usage.

    Important exclusions:
    - Path construction/inspection/setup helpers are allowlisted and ignored
      for lineage (for example: `Path(...)`, `os.path.*`, `.mkdir()`,
      `.exists()`, `.is_file()`, `.resolve()`, and glob helpers).
    - `importlib.*` calls are treated as import/module setup and ignored for
      lineage purposes.
    - `os.environ[...] = ...` assignments are treated as environment setup and
      ignored for lineage purposes.
    - `sys.path.*` calls (for example `sys.path.insert(...)`) are treated
      as import-path setup and ignored for lineage purposes.
    """
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
    tracked_path_keys = set(checker.tracked_paths)
    truly_untracked_paths = {
        path: lines
        for path, lines in checker.untracked_paths.items()
        if path not in tracked_path_keys
    }

    missing_lineage: list[str] = []

    # A script is fully tracked only if it explicitly opens and closes tracking.
    if not checker.has_track_call:
        missing_lineage.append("Missing ln.track() call in script.")

    if not checker.has_finish_call:
        missing_lineage.append("Missing ln.finish() call in script.")

    if truly_untracked_paths:
        # Point to concrete source lines where paths are used.
        for fname, lines in sorted(truly_untracked_paths.items()):
            lines_str = ", ".join(f"line {l}" for l in lines)
            missing_lineage.append(f"File not tracked in lamindb: {fname} ({lines_str})")


    is_fully_tracked = not missing_lineage
    return VerifyLineageResult(
        is_fully_tracked=is_fully_tracked,
        missing_lineage=tuple(missing_lineage),
    )

