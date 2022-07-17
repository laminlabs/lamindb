from pathlib import Path
from typing import Dict

import sqlmodel as sqm
from lamindb_schema.id import id_dobject

import lamindb as db
from lamindb._setup import load_settings

from .._logger import colors, logger
from ..admin.db import get_engine
from ..dev.file import store_file


def track_ingest(dobject_id):
    engine = get_engine()

    from nbproject import meta

    settings = load_settings()

    user_id = settings.user_id

    interface_id = meta.store.id

    with sqm.Session(engine) as session:
        track_do = db.model.track_do(
            type="ingest",
            user_id=user_id,
            interface_id=interface_id,
            dobject_id=dobject_id,
        )
        session.add(track_do)
        session.commit()
        session.refresh(track_do)

    settings._update_cloud_sqlite_file()

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

    def commit(self, integrity: bool = False, i_confirm_i_saved: bool = False):
        """Commit files for ingestion.

        Args:
            integrity: Check the integrity of the notebook.
            i_confirm_i_saved: Only relevant outside Jupyter Lab as a safeguard against
                losing the editor buffer content because of accidentally publishing.

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
        from nbproject import meta, publish
        from tabulate import tabulate  # type: ignore

        from lamindb.admin.db import insert

        settings = load_settings()
        logs = []

        for filepath, (dobject_id, dobject_v) in self.status.items():
            dobject_id = insert.dobject(
                filepath.stem,
                filepath.suffix,
                dobject_id=dobject_id,
                dobject_v=dobject_v,
            )

            dobject_storage_key = f"{dobject_id}-{dobject_v}{filepath.suffix}"
            store_file(filepath, dobject_storage_key)

            track_ingest(dobject_id)

            logs.append(
                [
                    f"{filepath.name} ({dobject_id}, {dobject_v})",
                    f"{meta.live.title!r} ({meta.store.id}, {meta.store.version})",
                    f"{settings.user_email} ({settings.user_id})",
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

        if not integrity:
            logger.warning(
                f"{colors.yellow('Consider using Jupyter Lab for ingesting data!')}\n"  # noqa
                "    Interactive notebook integrity checks are currently only supported on Jupyter Lab.\n"  # noqa
                "    Alternatively, manually save your notebook directly before calling `do.ingest.commit(..., integrity=True)`."  # noqa
            )
        publish(
            integrity=integrity,
            i_confirm_i_saved=i_confirm_i_saved,
            calling_statement="commit(",
        )


ingest = Ingest()
