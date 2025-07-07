from __future__ import annotations

import subprocess
import sys
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger

if TYPE_CHECKING:
    from lamindb.models import Run


def track_python_environment(run: Run) -> None:
    env_dir = ln_setup.settings.cache_dir / "environments" / f"run_{run.uid}"
    filepath = env_dir / "run_env_pip.txt"
    if not env_dir.exists():
        filepath.parent.mkdir(parents=True)
    # create a requirements.txt
    # we don't create a conda environment.yml mostly for its slowness
    try:
        with open(filepath, "w") as f:
            result = subprocess.run(
                [sys.executable, "-m", "pip", "freeze"],
                stdout=f,
            )
    except OSError as e:
        result = None
        logger.warning(f"could not run pip freeze with error {e}")
    if result is not None and result.returncode == 0:
        logger.info(f"tracked pip freeze > {str(filepath)}")
