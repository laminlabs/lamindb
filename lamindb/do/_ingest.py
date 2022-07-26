from pathlib import Path
from typing import Dict

import sqlmodel as sqm
from lndb_schema_core import id
from lndb_setup import settings

import lamindb as db

from .._logger import colors, logger
from ..dev.file import store_file


def track_ingest(dobject_id, dobject_v):
    from nbproject import meta

    user_id = settings.user.user_id

    jupynb_id = meta.store.id
    jupynb_v = meta.store.version

    with sqm.Session(settings.instance.db_engine()) as session:
        track_do = db.schema.core.track_do(
            type="ingest",
            user_id=user_id,
            jupynb_id=jupynb_id,
            jupynb_v=jupynb_v,
            dobject_id=dobject_id,
            dobject_v=dobject_v,
        )
        session.add(track_do)
        session.commit()
        session.refresh(track_do)

    settings.instance._update_cloud_sqlite_file()

    return track_do.id


class Ingest:
    """Ingest file."""

    def __init__(self) -> None:
        self._added: Dict = {}

    @property
    def status(self) -> dict:
        """Added files for ingestion."""
        return self._added

    def add(self, filepath, *, dobject_id=None, dobject_v="1"):
        """Add a file for ingestion.

        Args:
            filepath: The filepath.
            dobject_id: The dobject id.
            dobject_v: The dobject version.
        """
        filepath = Path(filepath)
        primary_key = (id.id_dobject() if dobject_id is None else dobject_id, dobject_v)
        self._added[filepath] = primary_key

    def commit(self, jupynb_v=None):
        """Commit files for ingestion.

        Args:
            jupynb_v: Notebook version to publish. Is automatically bumped if None.

        We primarily work with base62 IDs.

        ====== =========
        len_id n_entries
        ====== =========
        1      >6e+01
        2      >4e+03
        3      >2e+05
        4      >1e+07
        5      >9e+08
        6      >6e+10
        7      >4e+12
        8      >2e+14
        9      >1e+16
        12     >3e+21 (nbproject id)
        20     >7e+35 (~UUID)
        ====== =========
        """
        from nbproject import dev, meta, publish
        from tabulate import tabulate  # type: ignore

        from lamindb.dev.db import insert

        logs = []

        if meta.live.title is None:
            raise RuntimeError(
                "Can only ingest from notebook with title. Please set a title!"
            )

        jupynb_id = meta.store.id
        jupynb_v = dev.set_version(jupynb_v)  # version to be set in publish()
        jupynb_name = meta.live.title
        for filepath, (dobject_id, dobject_v) in self.status.items():
            dobject_id = insert.dobject(
                name=filepath.stem,
                file_suffix=filepath.suffix,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                jupynb_type="nbproject",
                dobject_id=dobject_id,
                dobject_v=dobject_v,
            )

            dobject_storage_key = f"{dobject_id}-{dobject_v}{filepath.suffix}"
            store_file(filepath, dobject_storage_key)

            track_ingest(dobject_id, dobject_v)

            logs.append(
                [
                    f"{filepath.name} ({dobject_id}, {dobject_v})",
                    f"{jupynb_name!r} ({jupynb_id}, {jupynb_v})",
                    f"{settings.user.user_email} ({settings.user.user_id})",
                ]
            )

        log_table = tabulate(
            logs,
            headers=[
                colors.green("Ingested File"),
                colors.blue("Notebook"),
                colors.purple("User"),
            ],
            tablefmt="pretty",
        )
        logger.success(f"{colors.bold('Ingested the following files')}:\n{log_table}")
        publish(calling_statement="commit(")


ingest = Ingest()
