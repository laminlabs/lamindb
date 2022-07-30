from pathlib import Path
from typing import Dict

import sqlmodel as sqm
from lndb_schema_core import id
from lndb_setup import settings

import lamindb as db

from .._logger import colors, logger
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_file_suffix
from ..meta import FeatureModel


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

    def add(
        self,
        dobject,
        *,
        feature_model=None,
        dobject_id=None,
        dobject_v="1",
    ):
        """Add a dobject or a file for ingestion.

        Args:
            dobject: An data object or filepath.
            feature_model: The data model that defines feature to annotate.
            dobject_id: The dobject id.
            dobject_v: The dobject version.
        """
        primary_key = (
            id.id_dobject() if dobject_id is None else dobject_id,
            dobject_v,
        )

        if isinstance(dobject, (Path, str)):
            # if a file is given, create a dobject
            # TODO: handle compressed files
            filepath = Path(dobject)
        else:
            # if in-memory object is given, return the cache path
            filekey = f"{primary_key[0]}-{primary_key[1]}{infer_file_suffix(dobject)}"
            filepath = settings.instance.storage.key_to_filepath(filekey)

        # skip feature annotation for store-only files
        if (filepath.suffix in [".fastq", ".fastqc", ".bam", ".sam", ".png"]) or (
            feature_model is None
        ):
            self._added[filepath] = primary_key
            return None

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
        self.integrity_flag = {
            "feature": fm.id_type,
            "n_mapped": n_mapped,
            "percent_mapped": round(n_mapped / n * 100, 1),
        }

        self._features[filepath] = (fm, df_curated)

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
            dobject_id = insert.dobject_from_jupynb(
                name=filepath.stem,
                file_suffix=filepath.suffix,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
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

            if self._features.get(filepath) is not None:
                fm, df_curated = self._features.get(filepath)
                fm.ingest(df_curated)

        # pretty logging info
        log_table = tabulate(
            logs,
            headers=[
                colors.green("Ingested file"),
                colors.blue("Notebook"),
                colors.purple("User"),
            ],
            tablefmt="pretty",
        )
        logger.success(f"{colors.bold('Ingested the following files')}:\n{log_table}")

        publish(calling_statement="commit(")


ingest = Ingest()
