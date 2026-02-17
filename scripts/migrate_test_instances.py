#!/usr/bin/env python3
"""Migrate all LaminDB instances used in lamindb tests.

For each instance: connect, run migrations, create storage snapshot.
Run from repo root with: python scripts/migrate_test_instances.py
"""

import subprocess
import sys

INSTANCES = [
    "laminlabs/lamin-site-assets",
    "laminlabs/lamin-dev",
    "laminlabs/lamindata",
    "laminlabs/cellxgene",
    "laminlabs/bionty-assets",
    "laminlabs/pertdata",
]


def run(cmd: str) -> None:
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        sys.exit(result.returncode)


def main() -> None:
    for instance in INSTANCES:
        print(f"=== Migrating {instance} ===")
        run(f"lamin connect {instance}")
        run("lamin migrate deploy")
        run("lamin io snapshot")
        print()

    print("Done. All test instances migrated and snapshotted.")


if __name__ == "__main__":
    main()
