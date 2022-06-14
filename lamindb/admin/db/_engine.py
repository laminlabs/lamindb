from pathlib import Path

import sqlmodel as sqm


def get_database_file() -> Path:
    from lamindb import settings

    database_filename = str(settings.storage_root.stem).lower()  # type: ignore
    return settings.storage_root / f"{database_filename}.lndb"  # type: ignore


def get_engine():
    database_file = get_database_file()
    return sqm.create_engine(f"sqlite:///{database_file}", future=True)
