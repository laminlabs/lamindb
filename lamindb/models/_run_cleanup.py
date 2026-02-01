"""Background cleanup of report/environment artifacts after Run bulk delete.

Runnable as: python -m lamindb.models._run_cleanup --instance owner/name --ids 1,2,3 [--run-uid UID]
"""

import argparse
import logging
from pathlib import Path

from lamin_utils import logger

import lamindb as ln


def main() -> None:
    parser = argparse.ArgumentParser(description="Clean up orphaned run artifacts.")
    parser.add_argument("--instance", required=True, help="Instance slug (owner/name).")
    parser.add_argument("--ids", required=True, help="Comma-separated artifact IDs.")
    parser.add_argument(
        "--run-uid",
        required=True,
        help="Run UID for log file name (run_cleanup_logs_{uid}.txt in cache dir).",
    )
    args = parser.parse_args()

    ln.connect(args.instance)

    file_handler = None
    log_path = (
        Path(ln.setup.settings.cache_dir) / f"run_cleanup_logs_{args.run_uid}.txt"
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_path, mode="a")
    logger.addHandler(file_handler)

    for aid_str in args.ids.split(","):
        aid = int(aid_str.strip())
        artifact = ln.Artifact.objects.filter(id=aid).first()
        if artifact is not None:
            assert artifact.kind == "__lamindb_run__", (
                f"artifact {artifact.uid} is not of __lamindb_run__ kind, aborting cleanup of artifacts {args.ids}"
            )
            try:
                artifact.delete(permanent=True)
                logger.important(f"deleted artifact {artifact.uid}")
            except Exception:
                logger.error(f"couldn't delete artifact {artifact.uid}")
                pass


if __name__ == "__main__":
    main()
