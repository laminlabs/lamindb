from pathlib import Path
from shutil import SameFileError
from typing import Dict

import sqlmodel as sqm
from lndb_schema_core import id
from lndb_setup import settings

import lamindb as db

from .._logger import colors, logger
from ..dev import storage_key_from_triple
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_file_suffix, write_to_file
from ._link import FeatureModel


def track_ingest(dobject_id, dobject_v):
    from nbproject import meta

    user_id = settings.user.id

    jupynb_id = meta.store.id
    jupynb_v = meta.store.version

    with sqm.Session(settings.instance.db_engine()) as session:
        usage = db.schema.core.usage(
            type="ingest",
            user_id=user_id,
            jupynb_id=jupynb_id,
            jupynb_v=jupynb_v,
            dobject_id=dobject_id,
            dobject_v=dobject_v,
        )
        session.add(usage)
        session.commit()
        session.refresh(usage)

    settings.instance._update_cloud_sqlite_file()

    return usage.id


class Ingest:
    """Ingest dobject."""

    def __init__(self) -> None:
        self._added: Dict = {}
        self._features: Dict = {}
        self._logs: Dict = {}

    @property
    def status(self) -> dict:
        """Added dobjects for ingestion."""
        return {k.as_posix(): v for k, v in self._added.items()}

    @property
    def logs(self) -> dict:
        """Logs of feature annotation."""
        return {k.as_posix(): v for k, v in self._logs.items()}

    def add(
        self,
        dobject,
        *,
        name=None,
        feature_model=None,
        dobject_id=None,
        dobject_v="1",
    ):
        """Stage a data object (in memory or file) for ingestion.

        Args:
            dobject: A data object in memory or filepath.
            name: A name. Required if passing in memory object.
            feature_model: Features to link during ingestion.
            dobject_id: The dobject id.
            dobject_v: The dobject version.
        """
        primary_key = (
            id.id_dobject() if dobject_id is None else dobject_id,
            dobject_v,
        )

        dmem = None
        if isinstance(dobject, (Path, str)):
            # if a file is given, create a dobject
            # TODO: handle compressed files
            filepath = Path(dobject)
        else:
            # if in-memory object is given, return the cache path
            dmem = dobject
            suffix = infer_file_suffix(dobject)
            if name is None:
                raise RuntimeError("Provide name if ingesting in memory data.")
            filepath = Path(f"{name}{suffix}")

        if feature_model is not None:
            # load file into memory
            if isinstance(dobject, (Path, str)):
                dmem = load_to_memory(dobject)
            else:
                dmem = dobject

            # curate features
            # looks for the id column, if none is found, will assume in the index
            try:
                df = getattr(dmem, "var")
            except AttributeError:
                df = dmem
            fm = FeatureModel(feature_model)
            df_curated = fm.curate(df)
            n = df_curated["__curated__"].count()
            n_mapped = df_curated["__curated__"].sum()
            self._logs[filepath] = {
                "feature": fm.id_type,
                "n_mapped": n_mapped,
                "percent_mapped": round(n_mapped / n * 100, 1),
                "unmapped": df_curated.index[~df_curated["__curated__"]],
            }

            self._features[filepath] = (fm, df_curated)

        self._added[filepath] = primary_key

        if not filepath.exists() and dmem is not None:
            write_to_file(dmem, filepath)

    def commit(self, jupynb_v=None):
        """Complete ingestion.

        Args:
            jupynb_v: Notebook version to publish. Is automatically set if `None`.
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
        for filepath, (dobject_id, dobject_v) in self._added.items():
            dobject_id = insert.dobject_from_jupynb(
                name=filepath.stem,
                file_suffix=filepath.suffix,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                dobject_v=dobject_v,
            )

            dobject_storage_key = storage_key_from_triple(
                dobject_id, dobject_v, filepath.suffix
            )
            try:
                store_file(filepath, dobject_storage_key)
            except SameFileError:
                pass

            track_ingest(dobject_id, dobject_v)

            logs.append(
                [
                    f"{filepath.name} ({dobject_id}, {dobject_v})",
                    f"{jupynb_name!r} ({jupynb_id}, {jupynb_v})",
                    f"{settings.user.handle} ({settings.user.id})",
                ]
            )

            if self._features.get(filepath) is not None:
                fm, df_curated = self._features.get(filepath)
                fm.ingest(dobject_id, df_curated)

        # pretty logging info
        log_table = tabulate(
            logs,
            headers=[
                colors.green("dobject"),
                colors.blue("jupynb"),
                colors.purple("user"),
            ],
            tablefmt="pretty",
        )
        logger.success(f"Ingested the following dobjects:\n{log_table}")

        publish(calling_statement="commit(")


ingest = Ingest()
