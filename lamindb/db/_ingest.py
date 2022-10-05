from pathlib import Path
from typing import Any, Dict, Optional, Union

from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id

from .._logger import colors, logger
from ..dev import get_name_suffix_from_filepath, track_usage
from ..dev.file import load_to_memory, store_file
from ..dev.object import infer_suffix, write_to_file
from ..schema import core, list_entities
from ._insert import insert
from ._link import link
from ._query import query


def store_insert_dobject(
    filepath: Path, dobject: core.dobject, dtransform: Union[core.jupynb, BfxRun]
):
    """Store and insert dobject."""
    dobject_storage_key = f"{dobject.id}{dobject.suffix}"
    size = store_file(filepath, dobject_storage_key)
    dobject.size = size

    if isinstance(dtransform, core.jupynb):
        insert.dobject_from_jupynb(dobject=dobject, jupynb=dtransform)  # type: ignore
        dtransform_log = dict(
            jupynb=f"{dtransform.name!r} ({dtransform.id}, {dtransform.v})"
        )
    elif isinstance(dtransform, BfxRun):
        insert.dobject_from_pipeline(  # type: ignore
            dobject=dobject, pipeline_run=dtransform
        )
        dtransform_log = dict(
            pipeline_run=f"{dtransform.run_name} ({dtransform.run_id})"
        )
    else:
        # TODO: implement non-jupynb insert
        raise NotImplementedError

    track_usage(dobject.id, usage_type="ingest")

    return dtransform_log


def insert_table_entries(
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


class Ingest:
    """Ingest dobjects.

    Ingest is an operation that stores, tracks and annotates dobjects.
    """

    def __init__(self) -> None:
        self._ingests: Dict = {}  # Ingest instances
        self._logs: list = []

        self._jupynb = None
        from nbproject import meta
        from nbproject.dev._jupyter_communicate import notebook_path

        if notebook_path() is not None:
            self._jupynb = core.jupynb(id=meta.store.id, name=meta.live.title)

    @property
    def userlog(self) -> dict:
        """User logging."""
        return dict(user=f"{settings.user.handle} ({settings.user.id})")

    @property
    def status(self) -> list:
        """Staged dobjects for ingestion."""
        dobjects = []
        for filepath_str, dobject_ in self._ingests.items():
            entry = dict(filepath=filepath_str, dobject_id=dobject_.dobject.id)
            dobjects.append(entry)
        return dobjects

    @property
    def ingests(self) -> dict:
        """Instances of ingest."""
        return self._ingests

    @property
    def jupynb(self):
        """Jupyter notebook."""
        return self._jupynb

    def _print_logging_table(
        self,
        message: str = "Ingested the following dobjects:",
    ):
        """Pretty print logs."""
        import pandas as pd
        from tabulate import tabulate  # type: ignore

        if len(self._logs) == 0:
            return

        log_table = tabulate(
            pd.DataFrame(self._logs).fillna(""),
            headers="keys",
            tablefmt="pretty",
            stralign="left",
        )

        logger.success(f"{message}\n{log_table}")
        """"""

    def add(self, data: Any, *, name: str = None, dobject_id: str = None):
        """Stage dobject for ingestion."""
        if isinstance(data, BfxRun):
            ingest_ = IngestPipelineRun(data)
        else:
            ingest_ = IngestDobject(  # type: ignore
                data, name=name, dobject_id=dobject_id
            )
        self._ingests[ingest_.filepath.as_posix()] = ingest_
        return ingest_

    def remove(self, filepath: Union[str, Path]):
        """Remove a dobject from the staged list."""
        filepath_str = filepath if isinstance(filepath, str) else filepath.as_posix()
        self._ingests.pop(filepath_str)

    def commit(self, jupynb_v=None):
        """Complete ingestion.

        Args:
            jupynb_v: Notebook version to publish. Is automatically set if `None`.
        """
        if self.jupynb is None:
            raise NotImplementedError
        else:
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
                    decide = input(
                        "   Do you still want to proceed with ingesting? (y/n) "
                    )

                if decide not in ("y", "Y", "yes", "Yes", "YES"):
                    logger.warning("Aborted!")
                    return "aborted"

            # version to be set in publish()
            self._jupynb.v = dev.set_version(jupynb_v)

            for filepath_str, ingest_ in self._ingests.items():
                # TODO: run the appropriate clean-up operations if any aspect
                # of the ingestion fails
                dtransformlog = ingest_._commit(dtransform=self.jupynb)

                self._logs.append({**ingest_.datalog, **dtransformlog, **self.userlog})

            self._print_logging_table()

            meta.store.version = self.jupynb.v
            meta.store.pypackage = meta.live.pypackage
            logger.info(
                f"Set notebook version to {colors.bold(meta.store.version)} & wrote"
                " pypackages."
            )
            meta.store.write(calling_statement="commit(")


class IngestDobject:
    """Ingest a dobject."""

    def __init__(self, data: Any, *, name: str = None, dobject_id: str = None) -> None:
        self._data = data  # input data object provided by user
        self._dmem = None  # in-memory object
        if isinstance(data, (Path, str)):
            # if a file is given, create a dobject
            # TODO: handle compressed files
            self._filepath = Path(data)
            name, suffix = get_name_suffix_from_filepath(self.filepath)
        else:
            # if in-memory object is given, return the cache path
            self._dmem = data
            suffix = infer_suffix(data)
            if name is None:
                raise RuntimeError("Provide name if ingesting in-memory data.")
            self._filepath = Path(f"{name}{suffix}")
            # write to file
            # TODO: should this raise an error/warning if filepath exists?
            if not self.filepath.exists():
                write_to_file(self.dmem, self.filepath)  # type: ignore

        # creates a dobject entry, but not inserted into the db yet
        dobject_id = id.id_dobject() if dobject_id is None else dobject_id
        self._dobject = core.dobject(id=dobject_id, name=name, suffix=suffix)

        # access to the feature model
        self._feature_model = None  # feature model

        # access to the link operations
        self._link = IngestLink(self)

    @property
    def data(self):
        """Data to be ingested."""
        return self._data

    @property
    def dobject(self):
        """An dobject entry to be inserted."""
        return self._dobject

    @property
    def link(self):
        """Link operations via ingest.

        See: `IngestLink`
        """
        return self._link

    @property
    def dmem(self) -> Any:
        """In-memory form of the dobject."""
        return self._dmem

    @property
    def filepath(self) -> Path:
        """Filepath of the dobject."""
        return self._filepath

    @property
    def feature_model(self):
        """Feature model."""
        return self._feature_model

    @property
    def datalog(self) -> dict:
        """Logging of a dobject entry.

        <filepath dobject_id>
        """
        return dict(dobject=f"{self.filepath.name} ({self.dobject.id})")

    def _commit(self, dtransform):
        """Store and insert dobject."""
        dtransform_log = store_insert_dobject(
            filepath=self.filepath, dobject=self.dobject, dtransform=dtransform
        )

        # insert features and link to dobject
        if self.feature_model is not None:
            self.feature_model["model"].ingest(
                self.dobject.id, self.feature_model["df_curated"]
            )

        return dtransform_log


class IngestPipelineRun:
    """Ingest pipeline runs.

    Data associated with pipeline runs are provided as "inputs" and "outputs".
    - "inputs" are directly associated with sample level metadata.
    - "outputs" are linked to inputs, through which also linked to metadata.
    """

    def __init__(self, run: BfxRun) -> None:
        self._run = run
        self._filepath = run.outdir
        self._ingests: Dict = {}  # IngestDobject instances
        for filepath in self.run.inputs + self.run.outputs:  # type: ignore
            self._ingests[filepath] = IngestDobject(filepath)

    @property
    def run(self):
        """Pipeline run."""
        return self._run

    @property
    def filepath(self) -> Path:
        """Filepath of the run."""
        return self._filepath

    @property
    def ingests(self) -> dict:
        """Added dobject instances for ingestion."""
        return self._ingests

    @property
    def datalog(self):
        """Logging of a dobject entry.

        <filepath dobject_id>
        """
        return dict(dobject=f"{self.filepath.name}")

    def _commit(self, dtransform=None):
        """Insert dobject and metadata entries."""
        dtransform = self.run

        # TODO: needs non-atomic erroring
        # register samples
        biometa_id = self._insert_biometa()

        # insert dobjects
        for _, ingest_ in self.ingests.items():
            dtransform_log = store_insert_dobject(
                filepath=ingest_.filepath,
                dobject=ingest_.dobject,
                dtransform=dtransform,
            )

        # insert metadata entries and link to dobjects
        self._link_biometa(biometa_id)
        self._link_pipeline_meta()

        # register entries in pipeline and pipeline_run tables
        insert_table_entries(
            table="pipeline",
            pk=self.run.pipeline_pk,
            name=self.run.pipeline_name,
            reference=self.run.pipeline_reference,
        )
        insert_table_entries(
            table="pipeline_run",
            pk=self.run.run_pk,
            fk={
                "pipeline_id": self.run.pipeline_pk["id"],
                "pipeline_v": self.run.pipeline_pk["v"],
            },
            name=self.run.outdir.as_posix(),
        )

        # register pipeline-specific pipeline and run tables
        # e.g. bfx_pipeline, bfx_run
        insert_table_entries(
            table=self.run.pipeline_table,
            pk=self.run.pipeline_pk,
        )
        insert_table_entries(
            table=self.run.run_table,
            pk=self.run.run_pk,
            fk=self.run.run_fk,
            dir=self.run.outdir.as_posix(),
        )

        return dtransform_log

    def _insert_biometa(self):
        """Register entries in bio metadata tables."""
        # insert biosample entries
        biosample_id = self.run.biosample_id
        if biosample_id is None:
            biosample_id = insert_table_entries(table="biosample", force=True)
        else:
            assert query.biosample(biosample_id=biosample_id).one()

        # insert techsample entries
        # TODO: insert needs to occur at the very end of commit
        # TODO: why force here?
        if "techsample" in list_entities():
            r1, r2 = self.run.fastq_path
            techsample_id = insert_table_entries(
                table="techsample",
                filepath_r1=r1.as_posix(),
                filepath_r2=r2.as_posix() if r2 is not None else r2,
                force=True,
            )
            # insert a biosample_techsample entry
            # TODO: needs to first check if the link entry exists
            insert_table_entries(
                table="biosample_techsample",
                biosample_id=biosample_id,
                techsample_id=techsample_id,
            )

        # insert a biometa entry
        # TODO: needs to check first before inserting
        biometa_id = insert_table_entries(
            table="biometa", biosample_id=biosample_id, force=True
        )

        return biometa_id

    def _link_biometa(self, biometa_id):
        """Link dobjects to a biometa entry."""
        for ingest_ in self.ingests.values():
            insert_table_entries(
                table="dobject_biometa",
                dobject_id=ingest_.dobject.id,
                biometa_id=biometa_id,
            )

    def _link_pipeline_meta(self):
        """Link dobjects to their pipeline-related metadata."""
        for filepath, ingest_ in self.ingests.items():
            # ingest pipeline-related metadata
            file_type = self.run.file_type.get(filepath)
            dir = filepath.parent.resolve().as_posix()
            pipeline_meta_id = insert_table_entries(
                table=self.run.meta_table, file_type=file_type, dir=dir
            )
            # link dobject to pipeline-related metadata
            insert_table_entries(
                table=f"dobject_{self.run.meta_table}",
                pk={
                    "dobject_id": ingest_.dobject.id,
                    f"{self.run.meta_table}_id": pipeline_meta_id,
                },
            )


class IngestLink:
    def __init__(self, dobject_: IngestDobject) -> None:
        self._dobject_ = dobject_

    def features(self, feature_model, *, featureset_name: str = None):
        # curate features
        # looks for the id column, if none is found, will assume in the index
        if self._dobject_.dmem is None:
            self._dobject_._dmem = load_to_memory(self._dobject_.data)
        try:
            df = getattr(self._dobject_.dmem, "var")  # for AnnData objects
            if callable(df):
                df = self._dobject_.dmem
        except AttributeError:
            df = self._dobject_.dmem

        self._dobject_._feature_model = link.feature_model(
            df=df, feature_model=feature_model, featureset_name=featureset_name
        )


ingest = Ingest()
