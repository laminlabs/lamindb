from pathlib import Path
from subprocess import run

import lndb
import pytest
from lamin_logger import logger


def pytest_sessionstart(session: pytest.Session):
    instance_dirs = [
        d for d in ["./docs/guide/mydata", "./mydata-test-db"] if Path(d).exists()
    ]
    for instance_dir in instance_dirs:
        cmd = f"rm -r {instance_dir}"
        try:
            run(cmd)
            logger.info(cmd)
        except Exception:
            logger.info(f"Could not delete {instance_dir}")
    cmd = "lamin init --storage mydata-test-db"
    logger.info(cmd)
    run(cmd, shell=True)
    try:
        lndb.load("testuser1/lamindb-ci", migrate=True)
        lndb.delete("lamindb-ci")
    except Exception:
        logger.info("Could not delete testuser1/lamindb-ci")


# given we moved away from nox virtual envs, we don't need the below anymore
# def pytest_sessionfinish(session: pytest.Session):
#     import session_info

#     session_info.show(
#         dependencies=True,
#         html=False,
#         excludes=[
#             "builtins",
#             "stdlib_list",
#             "importlib_metadata",
#             # Special module present if test coverage being calculated
#             # https://gitlab.com/joelostblom/session_info/-/issues/10
#             "$coverage",
#         ],
#     )
