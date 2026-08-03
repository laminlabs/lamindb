"""Background cleanup of report/environment artifacts after Run bulk delete.

Runnable as: python -m lamindb.models._run_cleanup --instance owner/name --ids 1,2,3 [--run-uid UID]
"""

import argparse
import logging
import traceback

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

    # Parent redirects stdout/stderr to the log file; also attach a FileHandler once
    # settings are available so lamin_utils logs are persisted even if streams change.
    print(f"run_cleanup start instance={args.instance} ids={args.ids}", flush=True)
    try:
        ln.connect(args.instance)
        log_path = ln.setup.settings.cache_dir / f"run_cleanup_logs_{args.run_uid}.txt"
        file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
        logger.addHandler(file_handler)
        logger.important(f"connected, cleaning up artifacts: {args.ids}")

        for aid_str in args.ids.split(","):
            aid = int(aid_str.strip())
            artifact = ln.Artifact.objects.filter(id=aid).first()
            if artifact is None:
                logger.warning(f"artifact id={aid} not found, skipping")
                continue
            assert artifact.kind == "__lamindb_run__", (
                f"artifact {artifact.uid} is not of __lamindb_run__ kind, aborting cleanup of artifacts {args.ids}"
            )
            try:
                artifact.delete(permanent=True)
                logger.important(f"deleted artifact {aid}")
            except Exception as e:
                logger.error(f"did not delete artifact {aid}: {e}")
        logger.important("run_cleanup finished")
    except Exception:
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
