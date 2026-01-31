"""Background cleanup of report/environment artifacts after Run bulk delete.

Runnable as: python -m lamindb.models._run_cleanup --instance owner/name --ids 1,2,3
"""

import argparse

from lamin_utils import logger

import lamindb as ln


def main() -> None:
    parser = argparse.ArgumentParser(description="Clean up orphaned run artifacts.")
    parser.add_argument("--instance", required=True, help="Instance slug (owner/name).")
    parser.add_argument("--ids", required=True, help="Comma-separated artifact IDs.")
    args = parser.parse_args()

    ln.connect(args.instance)

    for aid_str in args.ids.split(","):
        aid = int(aid_str.strip())
        artifact = ln.Artifact.objects.filter(id=aid).first()
        if artifact is None:
            continue
        assert artifact.kind == "__lamindb_run__", (
            f"artifact {artifact.uid} is not of __lamindb_run__ kind, aborting cleanup of artifacts {args.ids}"
        )
        try:
            artifact.delete(permanent=True)
        except Exception:
            logger.error(f"couldn't delete artifact {artifact.uid}")
            pass


if __name__ == "__main__":
    main()
