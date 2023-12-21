import subprocess

import lamindb_setup as ln_setup
from lamin_utils import logger
from lnschema_core.models import Run


def track_environment(run: Run):
    filepath = ln_setup.settings.storage.cache_dir / f"run_env_pip_{run.uid}"
    # create a requirements.txt
    result = subprocess.run(f"pip freeze > {str(filepath)}", shell=True)
    if result.returncode == 0:
        logger.info(f"tracked pip freeze > {str(filepath)}")
    # we don't create a conda environment.yml mostly for its slowness
    # result = subprocess.run(
    #     f"conda env export > {str(filepath)}", shell=True
    # )
    # if result.returncode == 0:
    #     with open(filepath) as f:
    #         content = f.read()
    #     with open(filepath, "w") as f:
    #         f.write("environment.yml\n{content}")
    #     logger.info("cached conda environment.yml")
