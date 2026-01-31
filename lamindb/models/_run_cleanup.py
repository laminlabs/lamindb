"""Background cleanup of report/environment artifacts after Run bulk delete.

Runnable as: python -m lamindb.models._run_cleanup --instance owner/name --ids 1,2,3
"""

from __future__ import annotations

import argparse
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description="Clean up orphaned run artifacts.")
    parser.add_argument("--instance", required=True, help="Instance slug (owner/name).")
    parser.add_argument("--ids", required=True, help="Comma-separated artifact IDs.")
    args = parser.parse_args()

    import lamindb as ln

    ln.connect(args.instance)

    from django.db.models import ProtectedError
    from django.db.utils import IntegrityError as DjangoIntegrityError

    from lamindb.errors import IntegrityError as LaminIntegrityError
    from lamindb.models import Artifact

    for aid_str in args.ids.split(","):
        aid = int(aid_str.strip())
        artifact = Artifact.filter(id=aid).first()
        if artifact is None:
            continue
        try:
            artifact.delete(permanent=True, storage=False)
        except (DjangoIntegrityError, LaminIntegrityError, ProtectedError):
            pass


if __name__ == "__main__":
    main()
    sys.exit(0)
