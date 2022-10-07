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


class init_ingest:
    @classmethod
    def jupynb(cls) -> Optional[core.jupynb]:
        """Check if currently inside a Jupyter Notebook."""
        from nbproject import meta
        from nbproject.dev import notebook_path

        if notebook_path() is not None:
            return core.jupynb(id=meta.store.id, name=meta.live.title)
        return None

    @classmethod
    def ingests(cls) -> dict:
        return {}

    @classmethod
    def logs(cls) -> list:
        return []

    @classmethod
    def userlog(cls) -> dict:
        return dict(user=f"{settings.user.handle} ({settings.user.id})")


_ingests = init_ingest.ingests()  # Ingest instances
_logs = init_ingest.logs()  # logging messages
jupynb = init_ingest.jupynb()
userlog = init_ingest.userlog()


class ingest:
    """Ingest operations.

    Ingest is an operation that stores, tracks and annotates dobjects.
    """

    @classmethod
    def status(cls) -> list:
        """Staged dobjects for ingestion."""
        dobjects = []
        for filepath_str, dobject_ in cls.list_ingests().items():
            entry = dict(filepath=filepath_str, dobject_id=dobject_.dobject.id)
            dobjects.append(entry)
        return dobjects

    @classmethod
    def list_ingests(cls) -> dict:
        """Instances of ingest."""
        return _ingests

    @classmethod
    def print_logging_table(
        cls,
        message: str = "Ingested the following dobjects:",
    ):
        """Pretty print logs."""
        import pandas as pd
        from tabulate import tabulate  # type: ignore

        if len(_logs) == 0:
            return

        log_table = tabulate(
            pd.DataFrame(_logs).fillna(""),
            headers="keys",
            tablefmt="pretty",
            stralign="left",
        )

        logger.success(f"{message}\n{log_table}")
        """"""

    @classmethod
    def reset(cls):
        global _ingests, _logs, jupynb, userlog
        _ingests = init_ingest.ingests()  # Ingest instances
        _logs = init_ingest.logs()  # logging messages
        jupynb = init_ingest.jupynb()
        userlog = init_ingest.userlog()

    @classmethod
    def add(cls, data: Any, *, name: str = None, dobject_id: str = None):
        """Stage dobject for ingestion."""
        ingest = Ingest(data, name=name, dobject_id=dobject_id)
        _ingests[ingest.filepath.as_posix()] = ingest
        return ingest

    @classmethod
    def remove(cls, filepath: Union[str, Path]):
        """Remove a dobject from the staged list."""
        filepath_str = filepath if isinstance(filepath, str) else filepath.as_posix()
        _ingests.pop(filepath_str)

    @classmethod
    def commit(cls, jupynb_v=None):
        """Complete ingestion.

        Args:
            jupynb_v: Notebook version to publish. Is automatically set if `None`.
        """
        if jupynb is None:
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
            jupynb.v = dev.set_version(jupynb_v)

            for filepath_str, ingest in cls.list_ingests().items():
                # TODO: run the appropriate clean-up operations if any aspect
                # of the ingestion fails
                ingest.commit()

                _logs.append({**ingest.datalog, **ingest.dtransformlog, **userlog})

            cls.print_logging_table()

            meta.store.version = jupynb.v
            meta.store.pypackage = meta.live.pypackage
            logger.info(
                f"Set notebook version to {colors.bold(meta.store.version)} & wrote"
                " pypackages."
            )
            meta.store.write(calling_statement="commit(")

        # reset ingest
        cls.reset()


class Ingest:
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

        # dtransform
        self._dtransform = None

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

    @property
    def dtransformlog(self) -> Optional[dict]:
        """Logging of the dtransform."""
        return self._dtransformlog

    @property
    def dtransform(self) -> Optional[core.dtransform]:
        """An dtransform entry to be inserted."""
        return self._dtransform

    @dtransform.setter
    def dtransform(self, value: core.dtransform):
        """Set value of dtransform."""
        self._dtransform = value

    def commit(self):
        """Store and insert dobject."""
        if self.dtransform is None:
            if jupynb is not None:
                self.link.jupynb(jupynb)
                self._dtransformlog = dict(
                    jupynb=f"{jupynb.name!r} ({jupynb.id}, {jupynb.v})"
                )
            else:
                raise RuntimeError("dtransform can't be None!")

        # store dobject
        dobject_storage_key = f"{self.dobject.id}{self.dobject.suffix}"
        size = store_file(self.filepath, dobject_storage_key)
        self._dobject.size = size  # size is only calculated when storing the file

        # if isinstance(self.dtransform, core.pipeline_run):
        #     self.link.pipeline_run(self.dtransform)
        #     dtransform_log = dict(
        #         jupynb=f"{self.dtransform.name!r} {self.dtransform.id}"
        #     )
        # else:
        #     raise NotImplementedError

        # insert all linked entries including dtransform
        for table_name, entry in self.link.linked_entries.items():
            getattr(insert, table_name)(**entry.dict())

        # insert dobject with storage_id and dtransform_id
        insert.dobject_from_dtransform(
            dobject=self.dobject, dtransform_id=self.dtransform.id
        )

        # insert features and link to dobject
        if self.feature_model is not None:
            self.feature_model["model"].ingest(
                self.dobject.id, self.feature_model["df_curated"]
            )

        track_usage(self.dobject.id, usage_type="ingest")


class IngestPipelineRun:
    """Ingest pipeline runs.

    Data associated with pipeline runs are provided as "inputs" and "outputs".
    - "inputs" are directly associated with sample level metadata.
    - "outputs" are linked to inputs, through which also linked to metadata.
    """

    def __init__(self, run: BfxRun) -> None:
        self._run = run
        self._filepath = run.outdir
        self._ingests: Dict = {}  # Ingest instances
        for filepath in self.run.inputs + self.run.outputs:  # type: ignore
            self._ingests[filepath] = Ingest(filepath)

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

    def commit(self):
        """Insert dobject and metadata entries."""
        # TODO: needs non-atomic erroring
        # register samples
        biometa_id = self._insert_biometa()

        # insert dobjects
        for _, ingest in self.ingests.items():
            dtransform_log = ingest.commit()

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
        for ingest in self.ingests.values():
            insert_table_entries(
                table="dobject_biometa",
                dobject_id=ingest.dobject.id,
                biometa_id=biometa_id,
            )

    def _link_pipeline_meta(self):
        """Link dobjects to their pipeline-related metadata."""
        for filepath, ingest in self.ingests.items():
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
                    "dobject_id": ingest.dobject.id,
                    f"{self.run.meta_table}_id": pipeline_meta_id,
                },
            )


class IngestLink:
    def __init__(self, dobject_: Ingest) -> None:
        self._dobject_ = dobject_
        self._entries: Dict = {}  # staged entries to be inserted

    @property
    def linked_entries(self) -> dict:
        return self._entries

    def features(self, feature_model, *, featureset_name: str = None):
        """Link dobject to features."""
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

    def pipeline_run(self, pipeline_run: core.pipeline_run) -> core.dtransform:
        """Link dobject to a pipeline run.

        - Create a dtransform entry from pipeline_run (no insert)
        - Link dobject to dtransform
        """
        dtransform = query.dtransform(  # type: ignore
            pipeline_run_id=pipeline_run.id
        ).one_or_none()
        if dtransform is None:
            dtransform = core.dtransform(pipeline_run_id=pipeline_run.id)

        self._entries["dtransform"] = dtransform
        self._dobject_.dtransform = dtransform

    def jupynb(self, jupynb: core.jupynb) -> core.dtransform:
        """Link dobject to a jupynb.

        - Create a dtransform entry from jupynb (no insert)
        - Link dobject to dtransform
        """
        result = query.jupynb(id=jupynb.id, v=jupynb.v).one_or_none()  # type: ignore
        if result is None:
            jupynb = core.jupynb(id=jupynb.id, v=jupynb.v, name=jupynb.name)
            self._entries["jupynb"] = jupynb
            # logger.info(
            #     f"Added notebook {jupynb.name!r} ({jupynb.id}, {jupynb.v}) by"
            #     f" user {settings.user.handle}."
            # )
            # dtransform entry
            dtransform = core.dtransform(jupynb_id=jupynb.id, jupynb_v=jupynb.v)
        else:
            dtransform = query.dtransform(  # type: ignore
                jupynb_id=jupynb.id, jupynb_v=jupynb.v
            ).one()

        self._entries["dtransform"] = dtransform
        self._dobject_.dtransform = dtransform
