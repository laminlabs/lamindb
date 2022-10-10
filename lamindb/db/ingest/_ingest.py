from pathlib import Path
from typing import Any, Dict, Optional, Union

from lndb_setup import settings
from lnschema_core import id

from ..._logger import logger
from ...dev import get_name_suffix_from_filepath, track_usage
from ...dev.file import load_to_memory, store_file
from ...dev.object import infer_suffix, write_to_file
from ...schema import core
from .._insert import insert
from .._link import link
from .._query import query


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


class Ingest:
    """Ingest data objects, initiated upon :class:`~lamindb.db.ingest.add`.

    Args:
        data: filepath or in-memory objects
        name: name of the data object, required of an in-memory object is passed
        dobject_id: id of the dobject

    Guide: :doc:`/db/guide/ingest`.
    """

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
        self._link = LinkIngest(self)

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

        Access point to :class:`~lamindb.db.ingest.LinkIngest`
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
        """Feature model used to ingest the features of dobject.

        See :class:`~lamindb.db.ingest.LinkFeatureModel`
        """
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
        """The dtransform entry linked to the dobject."""
        return self._dtransform

    @dtransform.setter
    def dtransform(self, value: core.dtransform):
        """Set value of dtransform."""
        self._dtransform = value

    def commit(self):
        """Store and insert dobject and its linked entries."""
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


class LinkIngest:
    """Link db entries to the dobject, accessible `Ingest().link`."""

    def __init__(self, ingest: Ingest) -> None:
        self._ingest = ingest
        # TODO: need to allow multiple entries from the same table
        self._entries: Dict = {}  # staged entries to be inserted

    @property
    def linked_entries(self) -> dict:
        """Linked db entries of the dobject."""
        return self._entries

    def features(self, feature_model, *, featureset_name: str = None):
        """Link dobject to features."""
        # curate features
        # looks for the id column, if none is found, will assume in the index
        if self._ingest.dmem is None:
            self._ingest._dmem = load_to_memory(self._ingest.data)
        try:
            df = getattr(self._ingest.dmem, "var")  # for AnnData objects
            if callable(df):
                df = self._ingest.dmem
        except AttributeError:
            df = self._ingest.dmem
        # insert feature entries
        # TODO: needs to be staged without inserting here
        self._ingest._feature_model = link.feature_model(
            df=df, feature_model=feature_model, featureset_name=featureset_name
        )

    def pipeline_run(self, pipeline_run: core.pipeline_run) -> core.dtransform:
        """Link dobject to a pipeline run.

        - Create a dtransform entry from pipeline_run (no insert)
        - Link dobject to dtransform
        """
        result = query.pipeline_run(id=pipeline_run.id).one_or_none()  # type: ignore
        if result is None:
            self._entries["pipeline_run"] = pipeline_run
            dtransform = core.dtransform(pipeline_run_id=pipeline_run.id)
        else:
            dtransform = query.dtransform(  # type: ignore
                pipeline_run_id=pipeline_run.id
            ).one()

        self._entries["dtransform"] = dtransform
        self._ingest.dtransform = dtransform

    def jupynb(self, jupynb: core.jupynb) -> core.dtransform:
        """Link dobject to a jupynb.

        - Create a dtransform entry from jupynb (no insert)
        - Link dobject to dtransform
        """
        result = query.jupynb(id=jupynb.id, v=jupynb.v).one_or_none()  # type: ignore
        if result is None:
            self._entries["jupynb"] = jupynb
            dtransform = core.dtransform(jupynb_id=jupynb.id, jupynb_v=jupynb.v)
        else:
            dtransform = query.dtransform(  # type: ignore
                jupynb_id=jupynb.id, jupynb_v=jupynb.v
            ).one()

        self._entries["dtransform"] = dtransform
        self._ingest.dtransform = dtransform


def add(data: Any, *, name: str = None, dobject_id: str = None) -> Ingest:
    """Stage dobject for ingestion.

    Args:
        data: filepath or in-memory objects
        name: name of the data object, required of an in-memory object is passed
        dobject_id: id of the dobject
    """
    ingest = Ingest(data, name=name, dobject_id=dobject_id)
    _ingests[ingest.filepath.as_posix()] = ingest
    return ingest


def remove(filepath: Union[str, Path]):
    """Remove a dobject from the staged list.

    Args:
        filepath: file path of the data object, one of `.status()`
    """
    filepath_str = filepath if isinstance(filepath, str) else filepath.as_posix()
    _ingests.pop(filepath_str)


def status() -> list:
    """List staged dobjects for ingestion."""
    dobjects = []
    for filepath_str, ingest in list_ingests().items():
        entry = dict(filepath=filepath_str, dobject_id=ingest.dobject.id)
        dobjects.append(entry)
    return dobjects


def reset() -> None:
    """Reset ingest, clear all staged data objects."""
    global _ingests, _logs, jupynb, userlog
    _ingests = init_ingest.ingests()  # Ingest instances
    _logs = init_ingest.logs()  # logging messages
    jupynb = init_ingest.jupynb()
    userlog = init_ingest.userlog()


def commit(jupynb_v: str = None, i_confirm_i_saved: bool = False):
    """Complete ingestion.

    Args:
        jupynb_v: Notebook version to publish. Is automatically set if `None`.
        i_confirm_i_saved: Only relevant outside Jupyter Lab as a safeguard against
            losing the editor buffer content because of accidentally publishing.
    """
    if jupynb is None:
        raise NotImplementedError
    else:
        from nbproject import dev
        from nbproject._publish import finalize_publish, run_checks_for_publish

        result = run_checks_for_publish(
            calling_statement="commit(", i_confirm_i_saved=i_confirm_i_saved
        )
        if result != "checks-passed":
            return result

        # version to be set in finalize_publish()
        jupynb.v = dev.set_version(jupynb_v)

        for filepath_str, ingest in list_ingests().items():
            # TODO: run the appropriate clean-up operations if any aspect
            # of the ingestion fails
            ingest.commit()
            _logs.append({**ingest.datalog, **ingest.dtransformlog, **userlog})

        logger.info(
            f"Added notebook {jupynb.name!r} ({jupynb.id}, {jupynb.v}) by"
            f" user {settings.user.handle}."
        )

        print_logging_table()

        finalize_publish(version=jupynb_v, calling_statement="commit(")

    # reset ingest
    reset()


def list_ingests() -> dict:
    """Ingest objects created via `.add`."""
    return _ingests


def print_logging_table(message: str = "Ingested the following dobjects:"):
    """Pretty print logging messages."""
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
