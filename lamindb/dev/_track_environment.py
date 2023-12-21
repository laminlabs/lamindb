from subprocess import run
from lnschema_core.models import Run
import lamindb_setup as ln_setup
from lamin_utils import logger


def track_environment(run: Run):
    filepath_no_suffix = ln_setup.settings.storage.cache_dir / "run_env_{run.uid}"
    # create a conda environment.yml
    result = run(f"conda env export > {str(filepath_no_suffix)}_environment.yml", shell=True)
    if result.returncode == 0:
        logger.info("tracking conda environment")
    else:
        # create a requirements.txt
        result = run(f"pip freeze > {str(filepath_no_suffix)}_requirements.txt", shell=True)
        if result.returncode == 0:
            logger.info("tracking pip requirements.txt")
