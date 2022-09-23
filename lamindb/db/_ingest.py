from pathlib import Path
from shutil import SameFileError
from typing import Dict, List

from lnbfx import BfxRun
from lnbfx.dev import get_bfx_files_from_dir, parse_bfx_file_type
from lndb_setup import settings
from lnschema_core import id

from .._logger import colors, logger
from ..dev import format_pipeline_logs, storage_key_from_triple, track_usage
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_file_suffix, write_to_file
from ._insert import insert
from ._link import link
from ._query import query


class Ingest:
    """Ingest dobjects and pipeline runs."""

    def __init__(self) -> None:
        self._ingest_object = IngestObject()
        self._ingest_bfx: List = []

    @property
    def status(self) -> dict:
        """Added dobjects for ingestion."""
        dobjects = {}  # type: ignore
        ingest_entities = [self._ingest_object] + self._ingest_bfx
        for ingest in ingest_entities:
            ingest_dobjects = {k.as_posix(): v for k, v in ingest.dobjects.items()}
            dobjects = {**dobjects, **ingest_dobjects}
        return dobjects

    def add(
        self,
        dobject,
        *,
        name: str = None,
        feature_model=None,
        featureset_name: str = None,
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
        if dobject.__class__ == BfxRun:
            ingest = IngestBfxRun()
            self._ingest_bfx.append(ingest)
        else:
            ingest = self._ingest_object  # type: ignore

        ingest.add(
            dobject,
            name,
            feature_model,
            featureset_name,
            dobject_id,
            dobject_v,
        )

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

        ingest_entities = [self._ingest_object] + self._ingest_bfx
        for ingest in ingest_entities:
            ingest.commit(jupynb_id, jupynb_v, jupynb_name)
            ingest.log()

        publish(calling_statement="commit(")


class IngestObject:
    """Ingest dobjects."""

    def __init__(self) -> None:
        self._dobjects: Dict = {}
        self._features: Dict = {}
        self._ingestion_logs: List = []

    @property
    def dobjects(self) -> Dict:
        """Added dobjects for ingestion."""
        return {k: v for k, v in self._dobjects.items()}

    def add(
        self,
        dobject,
        name,
        feature_model,
        featureset_name,
        dobject_id,
        dobject_v,
    ):
        """Stage dobject for ingestion."""
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
            self._features[filepath], self._feature_logs[filepath] = link.feature_model(
                df=df, feature_model=feature_model, featureset_name=featureset_name
            )

        self._dobjects[filepath] = primary_key

        if not filepath.exists() and dmem is not None:
            write_to_file(dmem, filepath)  # type: ignore

        return filepath

    def commit(self, jupynb_id, jupynb_v, jupynb_name):
        """Ingest staged dobjects."""
        for filepath, (dobject_id, dobject_v) in self._dobjects.items():
            dobject_id = ingest_dobject(
                filepath=filepath,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                dobject_v=dobject_v,
                pipeline_run=None,
            )

            self._ingestion_logs.append(
                [
                    f"{filepath.name} ({dobject_id}, {dobject_v})",
                    f"{jupynb_name!r} ({jupynb_id}, {jupynb_v})",
                    f"{settings.user.handle} ({settings.user.id})",
                ]
            )

            if self._features.get(filepath) is not None:
                fm, df_curated = self._features.get(filepath)
                fm.ingest(dobject_id, df_curated)

    def log(self):
        """Pretty print logs."""
        from tabulate import tabulate  # type: ignore

        if len(self._ingestion_logs) == 0:
            return

        log_table = tabulate(
            self._ingestion_logs,
            headers=[
                colors.green("dobject"),
                colors.blue("jupynb"),
                colors.purple("user"),
            ],
            tablefmt="pretty",
            stralign="left",
        )

        logger.success(f"Ingested the following dobjects:\n{log_table}")


class IngestBfxRun:
    """Ingest bfx runs."""

    def __init__(self) -> None:
        self._run = None
        self._outputs: Dict = {}
        self._inputs: Dict = {}
        self._ingestion_logs: List = []

    @property
    def dobjects(self) -> Dict:
        """Added dobjects for ingestion."""
        dobjects = {**self._inputs, **self._outputs}
        return {k: v for k, v in dobjects.items()}

    def add(self, run: BfxRun, *args):
        """Stage pipeline run entities for ingestion."""
        # stage bfx run
        self._run = run

        # stage inputs
        for fastq_path in self._run.fastq_path:  # type: ignore
            self._inputs[fastq_path] = (id.id_dobject(), "1")

        # stage outputs
        outputs = get_bfx_files_from_dir(self._run.outdir)  # type: ignore
        for output_path in outputs:
            self._outputs[Path(output_path)] = (id.id_dobject(), "1")

    def commit(self, jupynb_id, jupynb_v, jupynb_name):
        """Ingest staged pipeline run entities."""
        # register samples
        biometa_id = self._register_samples()

        # register run inputs and outputs as dobjects and link them to metadata
        self._register_dobjects(jupynb_id, jupynb_v, jupynb_name)
        self._link_dobjects_biometa(biometa_id)
        self._link_dobjects_bfxmeta()

        # register core pipeline and pipeline run
        self._register(
            "pipeline",
            pk={"id": self._run.pipeline_id, "v": self._run.pipeline_v},
            name=self._run.pipeline_name,
            reference=self._run.pipeline_reference,
        )
        self._register(
            "pipeline_run",
            pk={"id": self._run.run_id},
            pipeline_id=self._run.pipeline_id,
            pipeline_v=self._run.pipeline_v,
        )

        # register bfx pipeline and bfx pipeline run
        self._register(
            "bfx_pipeline", pk={"id": self._run.pipeline_id, "v": self._run.pipeline_v}
        )
        self._register(
            "bfx_run",
            pk={"id": self._run.run_id},
            pipeline_id=self._run.pipeline_id,
            pipeline_v=self._run.pipeline_v,
        )

    def log(self):
        """Pretty print logs."""
        from tabulate import tabulate  # type: ignore

        if len(self._ingestion_logs) == 0:
            return

        self._ingestion_logs = format_pipeline_logs(self._ingestion_logs)
        log_table = tabulate(
            self._ingestion_logs,
            headers=[
                colors.green("dobject"),
                colors.blue("pipeline run"),
                colors.purple("user"),
            ],
            tablefmt="pretty",
            stralign="left",
        )

        logger.success(
            f"Ingested the following dobjects from pipeline runs:\n{log_table}"
        )

    def _link_dobjects_biometa(self, biometa_id):
        """Link dobjects to a biometa entry."""
        for dobject_id, _ in self.dobjects.values():
            insert.dobject_biometa(dobject_id=dobject_id, biometa_id=biometa_id)

    def _link_dobjects_bfxmeta(self):
        """Register dobject's bfxmeta entry and link it to the dobject."""
        for filepath, (dobject_id, _) in self.dobjects.items():
            bfxmeta_id = self._register_bfxmeta(filepath)
            insert.dobject_bfxmeta(dobject_id=dobject_id, bfxmeta_id=bfxmeta_id)

    def _register(self, schema_module, pk: dict, **kwargs):
        """Generic registration helper function."""
        result = getattr(query, schema_module)(**pk).first()
        if result is None:
            getattr(insert, schema_module)(**pk, **kwargs)

    def _register_bfxmeta(self, dobject_filepath: Path):
        """Register bfxmeta entry."""
        file_type = parse_bfx_file_type(dobject_filepath, from_dir=True)
        dirpath = dobject_filepath.parent.resolve().as_posix()
        bfxmeta_id = insert.bfxmeta(file_type=file_type, dir=dirpath)  # type: ignore
        return bfxmeta_id

    def _register_dobjects(self, jupynb_id, jupynb_v, jupynb_name):
        """Register staged dobjects."""
        for filepath, (dobject_id, dobject_v) in self.dobjects.items():
            dobject_id = ingest_dobject(
                filepath=filepath,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                jupynb_name=jupynb_name,
                dobject_id=dobject_id,
                dobject_v=dobject_v,
                pipeline_run=self._run,
            )

            log_file_name = str(filepath.relative_to(self._run.run_dir))
            self._ingestion_logs.append(
                [
                    f"{log_file_name} ({dobject_id}, {dobject_v})",
                    f"{self._run.run_name} ({self._run.run_id})",
                    f"{settings.user.handle} ({settings.user.id})",
                ]
            )

    def _register_samples(self):
        """Register entries in the techsample, biosample, biometa tables."""
        sample_id = self._run.sample_id
        if sample_id is None:
            # insert techsample, biosample, and biometa entries
            if len(self._run.fastq_path) == 1:
                techsample_id = insert.techsample(filepath_r1=self._run.fastq_path)
            else:
                techsample_id = insert.techsample(
                    filepath_r1=self._run.fastq_path[0],
                    filepath_r2=self._run.fastq_path[1],
                )
            biosample_id = insert.biosample(techsample_id=techsample_id)
            biometa_id = insert.biometa(biosample_id=biosample_id)
        else:
            # check that biosample entry exists
            assert query.biosample(biosample_id=sample_id).one()
            # query for existing biometa
            result = query.biometa(
                where={"biosample": {"biosample_id": sample_id}}
            ).all()
            if len(result) == 0:
                # insert entry in biometa if it does not yet exist
                biometa_id = insert.biometa(biosample_id=sample_id)
            elif len(result) == 1:
                biometa_id = result[0].id
            else:
                raise AssertionError(
                    f"Multiple biometa entries associated with biosample {sample_id}."
                )

        return biometa_id


def ingest_dobject(
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


ingest = Ingest()
