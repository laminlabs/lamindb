from pathlib import Path
from subprocess import run

import lndb
from lamin_logger import logger


def pytest_sessionstart(session):
    instance_dirs = [
        d for d in ["./docs/guide/mydata", "./mydata-test-db"] if Path(d).exists()
    ]
    for instance_dir in instance_dirs:
        cmd = f"rm -r {instance_dir}"
        run(cmd)
        logger.info(cmd)
    cmd = "lamin init --storage mydata-test-db"
    logger.info(cmd)
    run(cmd, shell=True)
    try:
        lndb.load("testuser1/lamindb-ci", migrate=True)
        lndb.delete("lamindb-ci")
    except Exception:
        logger.info("Could not delete testuser1/lamindb-ci")
