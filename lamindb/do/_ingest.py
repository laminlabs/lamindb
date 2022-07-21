from pathlib import Path
from typing import Dict

import sqlmodel as sqm
from lamindb_schema.id import id_dobject

import lamindb as db
from lamindb._setup import (
    load_or_create_instance_settings,
    load_or_create_user_settings,
)

from .._logger import colors, logger
from ..admin.db import get_engine
from ..dev.file import store_file


def track_ingest(dobject_id):
    engine = get_engine()

    from nbproject import meta

    user_settings = load_or_create_user_settings()
    instance_settings = load_or_create_instance_settings()

    user_id = user_settings.user_id

    interface_id = meta.store.id

    with sqm.Session(engine) as session:
        track_do = db.schema.track_do(
            type="ingest",
            user_id=user_id,
            interface_id=interface_id,
            dobject_id=dobject_id,
        )
        session.add(track_do)
        session.commit()
        session.refresh(track_do)

    instance_settings._update_cloud_sqlite_file()

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
        primary_key = (id_dobject() if dobject_id is None else dobject_id, dobject_v)
        self._added[filepath] = primary_key

    def commit(self, interface_v=None):
        """Commit files for ingestion.

        Args:
            interface_v: Notebook version to publish. Is automatically bumped if None.

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

        from lamindb.admin.db import insert

        user_settings = load_or_create_user_settings()
        logs = []

        if meta.live.title is None:
            raise RuntimeError(
                "Can only ingest from notebook with title. Please set a title!"
            )

        interface_id = meta.store.id
        interface_v = dev.set_version(interface_v)  # version to be set in publish()
        interface_name = meta.live.title
        for filepath, (dobject_id, dobject_v) in self.status.items():
            dobject_id = insert.dobject(
                name=filepath.stem,
                file_suffix=filepath.suffix,
                interface_id=interface_id,
                interface_v=interface_v,
                interface_name=interface_name,
                interface_type="nbproject",
                dobject_id=dobject_id,
                dobject_v=dobject_v,
            )

            dobject_storage_key = f"{dobject_id}-{dobject_v}{filepath.suffix}"
            store_file(filepath, dobject_storage_key)

            track_ingest(dobject_id)

            logs.append(
                [
                    f"{filepath.name} ({dobject_id}, {dobject_v})",
                    f"{interface_name!r} ({interface_id}, {interface_v})",
                    f"{user_settings.user_email} ({user_settings.user_id})",
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
        logger.log(
            "INGEST", f"{colors.bold('Ingested the following files')}:\n{log_table}"
        )
        publish(calling_statement="commit(")


ingest = Ingest()
