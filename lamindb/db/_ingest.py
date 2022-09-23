from pathlib import Path
from shutil import SameFileError
from typing import Dict, List

from lnbfx import BfxRun
from lnbfx.dev import get_bfx_files_from_dir, parse_bfx_file_type
from lndb_setup import settings
from lnschema_core import id

from lamindb.dev._core import get_name_suffix_from_filepath

from .._logger import colors, logger
from ..dev import format_pipeline_logs, storage_key_from_triple, track_usage
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_suffix, write_to_file
from ._insert import insert
from ._link import link
from ._query import query


class Ingest:
    """Ingest dobjects and pipeline runs."""

    def __init__(self) -> None:
        self._ingest_object = IngestObject()  # dobjects
        self._ingest_bfx: List = []  # bfx pipeline runs

    @property
    def status(self) -> list:
        """Added dobjects for ingestion."""
        added_dobjects = {}  # type: ignore
        ingest_entities = self._ingest_bfx[:]
        ingest_entities.insert(0, self._ingest_object)
        for ingest in ingest_entities:
            added_dobjects = {**added_dobjects, **ingest.dobjects}

        added_list = []
        for k, v in added_dobjects.items():
            entry = dict(filepath=k.as_posix(), dobject_id=v[0], dobject_v=v[1])
            if k in self._ingest_object.feature_models.keys():
                entry["feature_model"] = self._ingest_object.feature_models[k]
            added_list.append(entry)

        return added_list

    @property
    def feature_models(self) -> dict:
        """Added feature models."""
        return self._ingest_object.feature_models

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
        self._feature_models: Dict = {}
        self._ingestion_logs: List = []

    @property
    def dobjects(self) -> Dict:
        """Added dobjects for ingestion."""
        return self._dobjects

    @property
    def feature_models(self) -> Dict:
        """Added feature models."""
        return self._feature_models

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

            # ingest with feature models
            curated = self._feature_models.get(filepath)
            if curated is not None:
                curated["model"].ingest(dobject_id, curated["df_curated"])

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
        return dobjects

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

            log_file_name = str(filepath.relative_to(self._run.outdir))
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
                techsample_id = insert.techsample(
                    filepath_r1=self._run.fastq_path[0].as_posix()
                )
            else:
                techsample_id = insert.techsample(
                    filepath_r1=self._run.fastq_path[0].as_posix(),
                    filepath_r2=self._run.fastq_path[1].as_posix(),
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
    name, suffix = get_name_suffix_from_filepath(filepath)

    if pipeline_run is None:
        dobject_id = insert.dobject_from_jupynb(
            name=name,
            suffix=suffix,
            jupynb_id=jupynb_id,
            jupynb_v=jupynb_v,
            jupynb_name=jupynb_name,
            dobject_id=dobject_id,
            dobject_v=dobject_v,
        )
    else:
        dobject_id = insert.dobject_from_pipeline(
            name=name,
            suffix=suffix,
            dobject_id=dobject_id,
            dobject_v=dobject_v,
            pipeline_run=pipeline_run,
        )

    dobject_storage_key = storage_key_from_triple(dobject_id, dobject_v, suffix)
    try:
        store_file(filepath, dobject_storage_key)
    except SameFileError:
        pass

    track_usage(dobject_id, dobject_v, usage_type="ingest")

    return dobject_id


ingest = Ingest()
