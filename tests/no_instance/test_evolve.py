import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "statement",
    [
        "from lamindb.setup import evolve",
        "import lamindb.setup.evolve as evolve",
        "evolve = ln.setup.evolve",
    ],
)
def test_evolve_imports_share_plan_declarations(tmp_path, statement):
    """The public module uses the backend's plan collector."""
    code = """
import sys
from uuid import UUID

def offline(event, args):
    if event in {"socket.connect", "socket.getaddrinfo"}:
        raise AssertionError("This import test must not access the network")

sys.addaudithook(offline)
from lamindb_setup import settings
from lamindb_setup.core import InstanceSettings
settings._instance_settings = InstanceSettings(
    id=UUID(int=0), owner="none", name="none", storage=None, modules="bionty,pertdb"
)
import lamindb as ln
assert not hasattr(ln, "evolve")
assert "lamindb_setup.evolve" not in sys.modules
"""
    code += (
        statement
        + """
import lamindb_setup.evolve as backend
assert evolve is backend
assert ln.setup.evolve is evolve

@evolve.create_plan
def measurements(ctx: backend.PlanContext, exports: None) -> None:
    backend.FeatureType("Measurements")

measurements.finalize()
assert [dep.obj.name for dep in measurements.own_dependencies] == ["Measurements"]
"""
    )
    env = {
        **os.environ,
        "LAMIN_SETTINGS_DIR": str(tmp_path),
        "LAMIN_ENV": "pytest",
        "PYTHONPATH": os.pathsep.join(
            [str(Path(__file__).resolve().parents[2]), os.environ.get("PYTHONPATH", "")]
        ),
    }
    env.pop("LAMIN_CURRENT_INSTANCE", None)
    env.pop("LAMIN_API_KEY", None)
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
