from pathlib import Path
from shutil import SameFileError
from typing import Dict

from lnbfx import BfxRun, get_bfx_files_from_dir
from lndb_setup import settings
from lnschema_core import id

from .._logger import colors, logger
from ..dev import storage_key_from_triple, track_usage
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_file_suffix, write_to_file
from ._insert import insert
from ._link import link
from ._query import query


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
        name: str = None,
        feature_model=None,
        featureset_name: str = None,
        pipeline_run: BfxRun = None,
        dobject_id=None,
        dobject_v="1",
    ):
        """Stage a data object (in memory or file) for ingestion.

        Args:
            dobject: A data object in memory or filepath.
            name: A name. Required if passing in memory object.
            feature_model: Features to link during ingestion.
            featureset_name: Name of the featureset to be linked.
            pipeline_run: The instance of pipeline run, e.g. BfxRun.
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
                if callable(df):
                    df = dmem
            except AttributeError:
                df = dmem
            self._features[filepath], self._logs[filepath] = link.feature_model(
                df=df, feature_model=feature_model, featureset_name=featureset_name
            )

        self._added[filepath] = primary_key

        # stage pipeline runs and run dobjects for ingestion
        if pipeline_run is not None:
            self._setup_ingestion_from_pipeline(filepath, pipeline_run)

        if not filepath.exists() and dmem is not None:
            write_to_file(dmem, filepath)  # type: ignore

    def commit(self, jupynb_v=None):
        """Complete ingestion.

        Args:
            jupynb_v: Notebook version to publish. Is automatically set if `None`.
        """
        from nbproject import dev, meta, publish

        if meta.live.title is None:
            raise RuntimeError(
                "Can only ingest from notebook with title. Please set a title!"
            )

        jupynb_id = meta.store.id
        jupynb_v = dev.set_version(jupynb_v)  # version to be set in publish()
        jupynb_name = meta.live.title

        # ingest pipeline entities (pipeline runs, pipelines, dobjects from pipelines)
        pipeline_logs = {}
        for run in set(self._pipeline_runs.values()):
            run_logs = []
            # ingest pipeline run and its pipeline
            self._ingest_pipeline_run(run)
            # ingest dobjects from the pipeline run
            dobject_paths = [
                path
                for path, pipeline_run in self._pipeline_runs.items()
                if pipeline_run is run
            ]
            for filepath in dobject_paths:
                dobject_id, dobject_v = self._added.get(filepath)
                dobject_id = self._ingest_dobject(
                    filepath=filepath,
                    jupynb_id=jupynb_id,
                    jupynb_v=jupynb_v,
                    jupynb_name=jupynb_name,
                    dobject_id=dobject_id,
                    dobject_v=dobject_v,
                    pipeline_run=run,
                )
                # link dobject to its bionformatics metadata (bfxmeta)
                run.link_dobject(dobject_id, filepath)

                log_file_name = str(filepath.relative_to(run.run_dir))
                run_logs.append(
                    [
                        f"{log_file_name} ({dobject_id}, {dobject_v})",
                        f"{jupynb_name!r} ({jupynb_id}, {jupynb_v})",
                        f"{settings.user.handle} ({settings.user.id})",
                    ]
                )

            pipeline_logs[(run.run_id, run.run_dir)] = run_logs

        # ingest jupynb entities (dobjects)
        jupynb_logs = []
        for filepath, (dobject_id, dobject_v) in self._added.items():
            # skip dobject if linked to a pipeline (already ingested)
            if self._pipeline_runs.get(filepath) is not None:
                continue
            dobject_id = self._ingest_dobject(
                filepath=filepath,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                dobject_v=dobject_v,
                pipeline_run=None,
            )

            jupynb_logs.append(
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
        self._log_from_pipeline(pipeline_logs)
        self._log_from_jupynb(jupynb_logs)

        publish(calling_statement="commit(")

    def _ingest_dobject(
        self,
        filepath,
        jupynb_id,
        jupynb_v,
        jupynb_name,
        dobject_id,
        dobject_v,
        pipeline_run,
    ):
        """Insert and store dobject."""
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

        dobject_storage_key = storage_key_from_triple(
            dobject_id, dobject_v, filepath.suffix
        )
        try:
            store_file(filepath, dobject_storage_key)
        except SameFileError:
            pass

        track_usage(dobject_id, dobject_v, "ingest")

        return dobject_id

    def _setup_ingestion_from_pipeline(self, input_path, pipeline_run):
        """Setup and stage pipeline run and run dobjects for ingestion."""
        # stage pipeline run for ingestion
        pipeline_run.db_engine = settings.instance.db_engine()
        self._pipeline_runs[input_path] = pipeline_run
        # parse pipeline run directory
        if input_path.is_dir():
            del self._added[input_path]  # do not ingest the directory
            del self._pipeline_runs[input_path]
            pipeline_run.run_dir = input_path
            dobjects_to_add = get_bfx_files_from_dir(input_path)
            for dobject in dobjects_to_add:
                self.add(dobject, pipeline_run=pipeline_run)

    def _ingest_pipeline_run(self, run):
        """Ingest staged pipeline runs and their pipelines."""
        from ._insert import insert

        # check if core pipeline exists, insert if not
        df_pipeline = getattr(query, "pipeline")(as_df=True).all()
        if (run.pipeline_id, run.pipeline_v) not in df_pipeline.index:
            insert.pipeline(
                id=run.pipeline_id,
                v=run.pipeline_v,
                name=run.pipeline_name,
                reference=run.pipeline_reference,
            )
        # check if core pipeline run exists, insert if not
        df_pipeline_run = getattr(query, "pipeline_run")(as_df=True).all()
        if run.run_id not in df_pipeline_run.index:
            insert.pipeline_run(
                id=run.run_id,
                name=run.run_name,
                pipeline_id=run.pipeline_id,
                pipeline_v=run.pipeline_v,
            )
        # check if bfx pipeline and run exist, insert if not
        run.check_and_ingest()

    def _log_from_pipeline(self, pipeline_logs: dict):
        """Pretty print logs from pipeline-related ingestion."""
        for (run_id, run_dir), logs in pipeline_logs.items():
            log_table = self._create_log_table(logs)
            logger.success(
                "Ingested the following dobjects from"
                f" pipeline run {run_id},"
                f" directory {run_dir}:"
                f"\n{log_table}"
            )

    def _log_from_jupynb(self, logs):
        """Pretty print logs from jupynb-related ingestion."""
        log_table = self._create_log_table(logs)
        logger.success(f"Ingested the following dobjects:\n{log_table}")

    def _create_log_table(self, logs):
        """Create log table for pretty printing."""
        from tabulate import tabulate  # type: ignore

        log_table = tabulate(
            logs,
            headers=[
                colors.green("dobject"),
                colors.blue("jupynb"),
                colors.purple("user"),
            ],
            tablefmt="pretty",
        )
        return log_table


ingest = Ingest()
