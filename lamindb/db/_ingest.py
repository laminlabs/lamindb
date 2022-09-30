from pathlib import Path
from shutil import SameFileError
from typing import Dict, Optional

from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id

from lamindb.dev._core import get_name_suffix_from_filepath

from .._logger import colors, logger
from ..dev import format_pipeline_logs, storage_key_from_triple, track_usage
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_suffix, write_to_file
from ..schema import list_entities
from ._insert import insert
from ._link import link
from ._query import query


class Ingest:
    """Ingest dobjects and pipeline runs."""

    def __init__(self) -> None:
        self._ingest_object = IngestObject()  # an ingest instance
        self._ingest_pipelines: list = []  # a list of ingest instances

    @property
    def status(self) -> list:
        """Added dobjects for ingestion."""
        added_dobjects = {}  # type: ignore
        ingest_entities = self._ingest_pipelines[:]
        ingest_entities.insert(0, self._ingest_object)
        for ingest in ingest_entities:
            added_dobjects = {**added_dobjects, **ingest.dobjects}

        added_list = []
        for k, v in added_dobjects.items():
            entry = dict(filepath=k.as_posix(), dobject_id=v)
            if k in self._ingest_object.feature_models.keys():
                entry["feature_model"] = self._ingest_object.feature_models[k]
            added_list.append(entry)

        return added_list

    def add(
        self,
        dobject,
        *,
        name: str = None,
        feature_model=None,
        featureset_name: str = None,
        dobject_id: str = None,
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
        if isinstance(dobject, BfxRun):
            ingest = IngestPipelineRun()
            self._ingest_pipelines.append(ingest)
        else:
            ingest = self._ingest_object  # type: ignore

        ingest.add(
            dobject,
            name,
            feature_model,
            featureset_name,
            dobject_id,
        )

    def commit(self, jupynb_v=None):
        """Complete ingestion.

        Args:
            jupynb_v: Notebook version to publish. Is automatically set if `None`.
        """
        from nbproject import dev, meta
        from nbproject.dev._check_last_cell import check_last_cell

        if meta.live.title is None:
            raise RuntimeError(
                "Can only ingest from notebook with title. Please set a title!"
            )
        nb = dev.read_notebook(meta._filepath)  # type: ignore
        if not check_last_cell(nb, calling_statement="commit("):
            raise RuntimeError(
                "Can only ingest from the last code cell of the notebook."
            )
        if not dev.check_consecutiveness(nb):
            if meta.env == "test":
                decide = "y"
            else:
                decide = input("   Do you still want to proceed with ingesting? (y/n) ")

            if decide not in ("y", "Y", "yes", "Yes", "YES"):
                logger.warning("Aborted!")
                return "aborted"

        jupynb_id = meta.store.id
        jupynb_v = dev.set_version(jupynb_v)  # version to be set in publish()
        jupynb_name = meta.live.title

        ingest_entities = self._ingest_pipelines[:]
        ingest_entities.insert(0, self._ingest_object)
        for ingest in ingest_entities:
            # TODO: run the appropriate clean-up operations if any aspect
            # of the ingestion fails
            ingest.commit(jupynb_id, jupynb_v, jupynb_name)
            ingest.log()

        meta.store.version = dev.set_version(jupynb_v)
        meta.store.pypackage = meta.live.pypackage
        logger.info(
            f"Set notebook version to {colors.bold(meta.store.version)} & wrote"
            " pypackages."
        )
        meta.store.write(calling_statement="commit(")


class IngestEntity:
    """Ingest any entity."""

    def __init__(self) -> None:
        self._dobjects: Dict = {}
        self._logs: list = []

    @property
    def dobjects(self) -> dict:
        """Added dobjects for ingestion."""
        return self._dobjects

    @property
    def logs(self):
        return self._logs

    def log(self, cols: tuple, message: str):
        """Pretty print logs."""
        from tabulate import tabulate  # type: ignore

        if len(self.logs) == 0:
            return

        log_table = tabulate(
            self.logs,
            headers=[
                colors.green(cols[0]),
                colors.blue(cols[1]),
                colors.purple(cols[2]),
            ],
            tablefmt="pretty",
            stralign="left",
        )

        logger.success(f"{message}\n{log_table}")

    def _ingest_dobject(
        self,
        filepath,
        jupynb_id,
        jupynb_v,
        jupynb_name,
        dobject_id,
        pipeline_run,
    ):
        """Helper function for dobject ingestion."""
        name, suffix = get_name_suffix_from_filepath(filepath)
        if dobject_id is None:
            dobject_id = id.dobject()

        dobject_storage_key = storage_key_from_triple(dobject_id, suffix)

        try:
            size = store_file(filepath, dobject_storage_key)
        except SameFileError:
            pass

        if pipeline_run is None:
            dobject_id = insert.dobject_from_jupynb(
                name=name,
                suffix=suffix,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                size=size,
            )
        else:
            dobject_id = insert.dobject_from_pipeline(
                name=name,
                suffix=suffix,
                dobject_id=dobject_id,
                pipeline_run=pipeline_run,
                size=size,
            )

        track_usage(dobject_id, usage_type="ingest")

        return dobject_id

    def _ingest_non_dobject(
        self,
        table: str,
        pk: Optional[dict] = None,
        fk: Optional[dict] = None,
        force: bool = False,
        **kwargs,
    ):
        """Generic helper function for non-dobject ingestion."""
        if pk is None:
            pk = {}
        if fk is None:
            fk = {}
        entry_id = getattr(insert, table)(**pk, **fk, **kwargs, force=force)
        if force is False and entry_id is None:
            entry_id = getattr(query, table)(**pk, **fk, **kwargs).one().id
        return entry_id


class IngestObject(IngestEntity):
    """Ingest dobjects."""

    def __init__(self) -> None:
        super().__init__()
        self._feature_models: Dict = {}

    @property
    def dobjects(self) -> dict:
        """Added dobjects for ingestion."""
        return self._dobjects

    @property
    def logs(self):
        return self._logs

    @property
    def feature_models(self) -> dict:
        """Added feature models."""
        return self._feature_models

    def add(
        self,
        dobject,
        name: str = None,
        feature_model=None,
        featureset_name: str = None,
        dobject_id: str = None,
    ):
        """Stage dobject for ingestion."""
        if dobject_id is None:
            dobject_id = id.id_dobject()

        dmem = None
        if isinstance(dobject, (Path, str)):
            # if a file is given, create a dobject
            # TODO: handle compressed files
            filepath = Path(dobject)
        else:
            # if in-memory object is given, return the cache path
            dmem = dobject
            suffix = infer_suffix(dobject)
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
            curated = link.feature_model(
                df=df, feature_model=feature_model, featureset_name=featureset_name
            )
            self._feature_models[filepath] = curated

        self._dobjects[filepath] = dobject_id

        if not filepath.exists() and dmem is not None:
            write_to_file(dmem, filepath)  # type: ignore

        return filepath

    def commit(self, jupynb_id, jupynb_v, jupynb_name):
        """Ingest staged dobjects."""
        for filepath, dobject_id in self.dobjects.items():
            dobject_id = self._ingest_dobject(
                filepath=filepath,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                pipeline_run=None,
            )

            self._logs.append(
                [
                    f"{filepath.name} ({dobject_id})",
                    f"{jupynb_name!r} ({jupynb_id}, {jupynb_v})",
                    f"{settings.user.handle} ({settings.user.id})",
                ]
            )

            # ingest with feature models
            curated = self._feature_models.get(filepath)
            if curated is not None:
                curated["model"].ingest(dobject_id, curated["df_curated"])

    def log(self):
        """Pretty print logs."""
        super().log(
            cols=("dobject", "jupynb", "user"),
            message="Ingested the following dobjects:",
        )


class IngestPipelineRun:
    """Ingest runs from any pipeline."""

    def __init__(self):
        super().__init__()
        self._run = None
        self._inputs: Dict = {}
        self._outputs: Dict = {}

    @property
    def dobjects(self) -> dict:
        """Added dobjects for ingestion."""
        return {**self._inputs, **self._outputs}

    @property
    def logs(self):
        return format_pipeline_logs(self._logs)

    def add(self, run, *args):
        # add run for ingestion
        self._run = run

        # stage inputs
        for path in self._run.inputs:  # type: ignore
            self._inputs[path] = id.id_dobject()

        # stage outputs
        for path in self._run.outputs:  # type: ignore
            self._outputs[Path(path)] = id.id_dobject()

    def commit(self, jupynb_id, jupynb_v, jupynb_name):
        # register samples
        biometa_id = self._ingest_samples()

        # register run inputs and outputs as dobjects and link them to metadata
        self._ingest_dobjects(jupynb_id, jupynb_v, jupynb_name)
        self._link_biometa(biometa_id)
        self._link_pipeline_meta()

        # register core pipeline and pipeline run
        self._ingest_non_dobject(
            table="pipeline",
            pk=self._run.pipeline_pk,
            name=self._run.pipeline_name,
            reference=self._run.pipeline_reference,
        )
        self._ingest_non_dobject(
            table="pipeline_run",
            pk=self._run.run_pk,
            fk={
                "pipeline_id": self._run.pipeline_pk["id"],
                "pipeline_v": self._run.pipeline_pk["v"],
            },
            name=self._run.outdir.as_posix(),
        )

        # register pipeline-specific pipeline and run
        self._ingest_non_dobject(
            table=self._run.pipeline_table,
            pk=self._run.pipeline_pk,
        )
        self._ingest_non_dobject(
            table=self._run.run_table,
            pk=self._run.run_pk,
            fk=self._run.run_fk,
            dir=self._run.outdir.as_posix(),
        )

    def log(self):
        super().log(
            cols=("dobject", "pipeline run", "user"),
            message="Ingested the following dobjects from pipeline runs:",
        )

    def _ingest_dobjects(self, jupynb_id, jupynb_v, jupynb_name):
        """Register staged dobjects."""
        for filepath, dobject_id in self.dobjects.items():
            dobject_id = self._ingest_dobject(
                filepath=filepath,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                pipeline_run=self._run,
            )

            try:
                log_file_name = str(filepath.relative_to(self._run.outdir))
            except ValueError:
                log_file_name = filepath.name
            self._logs.append(
                [
                    f"{log_file_name} ({dobject_id})",
                    f"{self._run.run_name} ({self._run.run_id})",
                    f"{settings.user.handle} ({settings.user.id})",
                ]
            )

    def _ingest_samples(self):
        """Register entries in bio metadata tables."""
        # ingest techsample
        if "techsample" in list_entities():
            if len(self._run.fastq_path) == 1:
                techsample_id = self._ingest_non_dobject(
                    table="techsample",
                    filepath_r1=self._run.fastq_path[0].as_posix(),
                    force=True,
                )
            else:
                techsample_id = self._ingest_non_dobject(
                    table="techsample",
                    filepath_r1=self._run.fastq_path[0].as_posix(),
                    filepath_r2=self._run.fastq_path[1].as_posix(),
                    force=True,
                )

        # ingest biosample
        biosample_id = self._run.biosample_id
        if biosample_id is None:
            biosample_id = self._ingest_non_dobject(table="biosample", force=True)
        else:
            assert query.biosample(biosample_id=biosample_id).one()

        # ingest biosample_techsample
        self._ingest_non_dobject(
            table="biosample_techsample",
            biosample_id=biosample_id,
            techsample_id=techsample_id,
        )

        # ingest biometa
        biometa_id = self._ingest_non_dobject(
            table="biometa", biosample_id=biosample_id, force=True
        )

        return biometa_id

    def _link_biometa(self, biometa_id):
        """Link dobjects to a biometa entry."""
        for dobject_id in self.dobjects.values():
            self._ingest_non_dobject(
                table="dobject_biometa",
                dobject_id=dobject_id,
                biometa_id=biometa_id,
            )

    def _link_pipeline_meta(self):
        """Link dobjects to their pipeline-related metadata."""
        for filepath, dobject_id in self.dobjects.items():
            # ingest pipeline-related metadata
            file_type = self._run.file_type.get(filepath)
            dir = filepath.parent.resolve().as_posix()
            pipeline_meta_id = self._ingest_non_dobject(
                table=self._run.meta_table, file_type=file_type, dir=dir
            )
            # link dobject to pipeline-related metadata
            self._ingest_non_dobject(
                table=f"dobject_{self._run.meta_table}",
                pk={
                    "dobject_id": dobject_id,
                    f"{self._run.meta_table}_id": pipeline_meta_id,
                },
            )


ingest = Ingest()
