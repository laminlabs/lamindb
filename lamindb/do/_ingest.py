from pathlib import Path
from shutil import SameFileError
from typing import Dict

from lnbfx import get_bfx_files_from_folder
from lndb_setup import settings
from lnschema_core import id

from .._logger import colors, logger
from ..dev import storage_key_from_triple, track_usage
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_file_suffix, write_to_file
from ._link import link


class Ingest:
    """Ingest dobject."""

    def __init__(self) -> None:
        self._added: Dict = {}
        self._features: Dict = {}
        self._pipeline_runs: Dict = {}
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
        pipeline_run=None,
        dobject_id=None,
        dobject_v="1",
    ):
        """Stage a data object (in memory or file) for ingestion.

        Args:
            dobject: A data object in memory or filepath.
            name: A name. Required if passing in memory object.
            feature_model: Features to link during ingestion.
            pipeline_run: The instance of pipeline run, e.g. BfxRun
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

        # ingest features via feature models
        if feature_model is not None:
            # load file into memory
            if isinstance(dobject, (Path, str)):
                dmem = load_to_memory(dobject)
            else:
                dmem = dobject
            # curate features
            # looks for the id column, if none is found, will assume in the index
            try:
                df = getattr(dmem, "var")  # for AnnData objects
            except AttributeError:
                df = dmem
            self._features[filepath], self._logs[filepath] = link.feature_model(
                df=df, feature_model=feature_model
            )

        self._added[filepath] = primary_key

        # pipeline run
        if pipeline_run is not None:
            if Path(dobject).is_dir():
                del self._added[filepath]
                pipeline_run.bfx_pipeline_run_folder = Path(dobject)
                dobjects_to_add = get_bfx_files_from_folder(dobject)
                for dobject in dobjects_to_add:
                    self.add(dobject, pipeline_run=pipeline_run)
            pipeline_run.db_engine = settings.instance.db_engine()
            self._pipeline_runs[filepath] = pipeline_run

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
            pipeline_run = self._pipeline_runs.get(filepath)
            dobject_id = insert.dobject_from_jupynb(
                name=filepath.stem,
                file_suffix=filepath.suffix,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                dobject_v=dobject_v,
                pipeline_run=pipeline_run,
            )

            if pipeline_run is not None:
                pipeline_run.link_dobject_to_bfxmeta(dobject_id, filepath)

            dobject_storage_key = storage_key_from_triple(
                dobject_id, dobject_v, filepath.suffix
            )
            try:
                store_file(filepath, dobject_storage_key)
            except SameFileError:
                pass

            track_usage(dobject_id, dobject_v, "ingest")

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
