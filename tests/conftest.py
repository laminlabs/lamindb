from pathlib import Path
from subprocess import run

import lndb
from lamin_logger import logger


def pytest_sessionstart(session):
    instance_dirs = [
        d for d in ["./docs/guide/mydata", "./mydata-test-db"] if Path(d).exists()
    ]
    for instance_dir in instance_dirs:
        clean_instance = f"rm -r {instance_dir}"
        logger.info(clean_instance)
        run(*clean_instance.split(" "))
    cmd = "lamin init --storage mydata-test-db"
    run(cmd, shell=True)
    logger.info(cmd)
    try:
        lndb.load("testuser1/lamindb-ci", migrate=True)
        lndb.delete("lamindb-ci")
    except Exception:
        logger.info("Could not delete testuser1/lamindb-ci")
