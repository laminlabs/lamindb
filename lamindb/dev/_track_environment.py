import subprocess

import lamindb_setup as ln_setup
from lamin_utils import logger
from lnschema_core.models import Run


def track_environment(run: Run) -> None:
    filepath = ln_setup.settings.storage.cache_dir / f"run_env_pip_{run.uid}"
    # create a requirements.txt
    # we don't create a conda environment.yml mostly for its slowness
    result = subprocess.run(f"pip freeze > {str(filepath)}", shell=True)
    if result.returncode == 0:
        logger.info(f"tracked pip freeze > {str(filepath)}")
