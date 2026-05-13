import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_NAMES = ("anndata", "h5py", "pyarrow")
LIGHT_IMPORTS = {name: False for name in MODULE_NAMES}


PROBE_CASES = [
    (
        "storage package constants stay light",
        "import lamindb.core.storage as storage\n_ = storage.VALID_SUFFIXES\n_ = storage.delete_storage\n_ = storage.infer_filesystem",
        LIGHT_IMPORTS,
        (),
    ),
    (
        "storage object helpers stay light",
        "import lamindb.core.storage as storage\n_ = storage.infer_suffix\n_ = storage.write_to_disk",
        LIGHT_IMPORTS,
        (),
    ),
    (
        "loaders basic helpers stay light",
        "import lamindb.core.loaders as loaders\n_ = loaders.load_json\n_ = loaders.load_txt\n_ = loaders.load_html",
        LIGHT_IMPORTS,
        (),
    ),
    (
        "loaders tabular helpers stay light",
        "import lamindb.core.loaders as loaders\n_ = loaders.load_csv\n_ = loaders.load_parquet\n_ = loaders.load_tsv",
        LIGHT_IMPORTS,
        (),
    ),
    (
        "loaders optional-format helpers stay light",
        "import lamindb.core.loaders as loaders\n_ = loaders.load_h5ad\n_ = loaders.load_h5mu\n_ = loaders.load_zarr",
        LIGHT_IMPORTS,
        (),
    ),
    (
        "backed_access symbols stay light",
        "from lamindb.core.storage._backed_access import BackedAccessor, backed_access, _open_dataframe\n_ = BackedAccessor\n_ = backed_access\n_ = _open_dataframe",
        LIGHT_IMPORTS,
        (),
    ),
    (
        "objects module import stays light",
        "from lamindb.core.storage.objects import infer_suffix, write_to_disk\n_ = infer_suffix\n_ = write_to_disk",
        LIGHT_IMPORTS,
        (),
    ),
    (
        "backed_access pyarrow dataframe path stays anndata-free",
        "from upath import UPath\nimport pyarrow as pa\nimport pyarrow.parquet as pq\nfrom lamindb.core.storage._backed_access import backed_access\npath = UPath('test_import_side_effects.parquet')\npq.write_table(pa.table({'col': [1]}), path.as_posix())\ntry:\n    _ = backed_access(path, engine='pyarrow')\nfinally:\n    if path.exists():\n        path.unlink()",
        {"anndata": False, "h5py": False, "pyarrow": True},
        ("pyarrow",),
    ),
    (
        "backed_access polars dataframe path stays light",
        "from upath import UPath\nfrom lamindb.core.storage._backed_access import backed_access\npath = UPath('test_import_side_effects.csv')\nwith path.open('w') as f:\n    _ = f.write('col\\n1\\n')\ntry:\n    _ = backed_access(path, engine='polars')\nfinally:\n    if path.exists():\n        path.unlink()",
        LIGHT_IMPORTS,
        ("polars",),
    ),
]


def _probe_modules_loaded(code: str) -> dict[str, bool]:
    env = os.environ.copy()
    pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not pythonpath
        else os.pathsep.join([str(REPO_ROOT), pythonpath])
    )
    probe_lines = [
        "import json",
        "import sys",
        "",
        f"module_names = {MODULE_NAMES!r}",
        "result = {name: (name in sys.modules) for name in module_names}",
        code,
        'result.update({f"{name}_after": (name in sys.modules) for name in module_names})',
        "print(json.dumps(result))",
    ]
    probe = "\n".join(probe_lines)
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        check=True,
        capture_output=True,
        cwd=REPO_ROOT,
        env=env,
        text=True,
    )
    stdout_lines = [line for line in completed.stdout.splitlines() if line.strip()]
    return json.loads(stdout_lines[-1])


def _assert_modules(
    result: dict[str, bool], expected_after: dict[str, bool], label: str
):
    for module_name in MODULE_NAMES:
        assert result[module_name] is False, (
            f"{label}: {module_name} loaded before probe"
        )
        assert result[f"{module_name}_after"] is expected_after[module_name], (
            f"{label}: unexpected {module_name} import state"
        )


@pytest.mark.parametrize(
    ("label", "code", "expected_after", "required_modules"),
    PROBE_CASES,
)
def test_storage_import_side_effects(
    label: str,
    code: str,
    expected_after: dict[str, bool],
    required_modules: tuple[str, ...],
):
    missing_modules = [
        module_name
        for module_name in required_modules
        if importlib.util.find_spec(module_name) is None
    ]
    if missing_modules:
        pytest.skip(f"missing optional dependency: {', '.join(missing_modules)}")

    result = _probe_modules_loaded(code)
    _assert_modules(result, expected_after, label)
