from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import sqlmodel as sqm
from lndb_setup import settings
from nbproject import dev, meta

from ..._logger import logger
from ...dev import get_name_suffix_from_filepath, track_usage
from ...dev.file import store_file
from ...dev.object import infer_suffix, write_to_file
from ...schema import core
from .._insert import insert
from .._query import query
from ._lingstaged import LinkStaged


class Staged:
    """Staged data objects, initiated upon :class:`~lamindb.db.ingest.add`.

    Args:
        data: filepath or in-memory objects
        dtransform: The data transformation that links to the data source of the data object.
        name: name of the data object, required of an in-memory object is passed
        dobject_id: id of the dobject

    Guide: :doc:`/db/guide/ingest`.
    """

    def __init__(
        self,
        data: Any,
        *,
        dtransform: core.dtransform,
        name: str = None,
        dobject_id: str = None,
    ) -> None:
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
        self._dobject = core.dobject(name=name, suffix=suffix)
        self._dobject.id = dobject_id if dobject_id is not None else self.dobject.id

        # access to the feature model
        self._feature_model = None  # feature model

        # access to the link operations
        self._link = LinkStaged(self)

        # dtransform
        self._dtransform = dtransform

    @property
    def data(self) -> Any:
        """Data provided by the user upon init."""
        return self._data

    @property
    def dobject(self) -> core.dobject:
        """An dobject entry to be inserted."""
        return self._dobject

    @property
    def link(self) -> LinkStaged:
        """Link operations via ingest."""
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
    def dtransform(self) -> core.dtransform:
        """The dtransform entry linked to the dobject."""
        return self._dtransform

    def cleanup(self) -> None:
        """Clean up all linked entries."""
        if "dtransform" in self._link._entries:
            self._dtransform = None
        self._link._entries = {}

    def _commit(self) -> None:
        """Store and insert dobject and its linked entries."""
        # store dobject
        dobject_storage_key = f"{self.dobject.id}{self.dobject.suffix}"
        size = store_file(self.filepath, dobject_storage_key)
        self._dobject.size = size  # size is only calculated when storing the file

        # insert all linked entries including dtransform
        for table_name, entries in self.link.linked_entries.items():
            insert.from_list(table_name=table_name, entries=entries)  # type:ignore

        # insert dobject with storage_id and dtransform_id
        insert.dobject_from_dtransform(  # type:ignore
            dobject=self.dobject, dtransform_id=self.dtransform.id  # type:ignore
        )

        # insert features and link to dobject
        if self.feature_model is not None:
            self.feature_model["model"].ingest(
                self.dobject.id, self.feature_model["df_curated"]
            )

        track_usage(self.dobject.id, usage_type="ingest")


class Ingest:
    """Ingest data.

    Ingest is an operation that stores and annotates data objects.

    Guide: :doc:`/db/guide/ingest`.

    Args:
        dsource: A data source. If `None` assumes a Jupyter Notebook. In the
            default schema `jupynb` and `pipeline_run` are the two allowed data
            sources.

    Ingest helper classes:

    .. autosummary::
        :toctree: .

        Staged
        LinkStaged
        LinkFeatureModel
    """

    def _init_dtransform(self, dsource):
        if isinstance(dsource, core.pipeline_run):
            dtransform = query.dtransform(  # type: ignore
                pipeline_run_id=dsource.id
            ).one_or_none()
            if dtransform is None:
                dtransform = core.dtransform(pipeline_run_id=dsource.id)
            log = dict(pipeline_run=f"{dsource.name!r} ({dsource.id})")
        elif isinstance(dsource, core.jupynb):
            dtransform = query.dtransform(  # type: ignore
                jupynb_id=dsource.id, jupynb_v=dsource.v
            ).one_or_none()
            if dtransform is None:
                dtransform = core.dtransform(jupynb_id=dsource.id, jupynb_v=dsource.v)
            log = dict(jupynb=f"{dsource.name!r} ({dsource.id}, {dsource.v})")
        return dtransform, log

    def __init__(self, dsource: Union[core.jupynb, core.pipeline_run, None] = None):
        if dsource is None:
            if dev.notebook_path() is not None:
                dsource = core.jupynb(id=meta.store.id, name=meta.live.title)
            else:
                raise RuntimeError("Please provide a data source.")
        self._dsource = dsource  # data source (pipeline_run or jupynb)
        self._dtransform, self._dtransformlog = self._init_dtransform(dsource)
        self._staged: Dict = {}  # staged dobjects
        self._logs: List = []  # logging messages
        self._userlog = dict(user=f"{settings.user.handle} ({settings.user.id})")
        self._committed = False

    def add(self, data: Any, *, name: str = None, dobject_id: str = None) -> Staged:
        """Stage dobject for ingestion.

        Args:
            data: filepath or in-memory objects
            name: name of the data object, required of an in-memory object is passed
            dobject_id: id of the dobject
        """
        ingest = Staged(
            data, dtransform=self._dtransform, name=name, dobject_id=dobject_id
        )
        self._staged[ingest.filepath.as_posix()] = ingest
        return ingest

    def remove(self, filepath: Union[str, Path]) -> None:
        """Remove a dobject from the staged list.

        Args:
            filepath: file path of the data object, one of `.status()`
        """
        filepath_str = filepath if isinstance(filepath, str) else filepath.as_posix()
        self._staged.pop(filepath_str)

    def status(self) -> list:
        """List staged dobjects for ingestion."""
        dobjects = []
        for filepath_str, ingest in self._staged.items():
            entry = dict(filepath=filepath_str, dobject_id=ingest.dobject.id)
            dobjects.append(entry)
        return dobjects

    def commit(self, jupynb_v: str = None, i_confirm_i_saved: bool = False) -> None:
        """Commit data to object storage and database.

        Args:
            jupynb_v: Notebook version to publish. Is automatically set if `None`.
            i_confirm_i_saved: Only relevant outside Jupyter Lab as a safeguard against
                losing the editor buffer content because of accidentally publishing.
        """
        if self._committed:
            logger.error("Already committed")
            return None

        if isinstance(self._dsource, core.jupynb):
            from nbproject._publish import finalize_publish, run_checks_for_publish

            result = run_checks_for_publish(
                calling_statement="commit(", i_confirm_i_saved=i_confirm_i_saved
            )
            if result != "checks-passed":
                return result

            # version to be set in finalize_publish()
            self._dsource.v = dev.set_version(jupynb_v)

        # insert dsource and dtransform
        with sqm.Session(settings.instance.db_engine()) as session:
            session.add(self._dsource)
            session.commit()  # need to commit here for foreign key integrity
            session.add(self._dtransform)
            session.commit()

        for filepath_str, ingest in self._staged.items():
            # TODO: run the appropriate clean-up operations if any aspect
            # of the ingestion fails
            ingest._commit()
            self._logs.append(
                {**ingest.datalog, **self._dtransformlog, **self._userlog}
            )

        if isinstance(self._dsource, core.jupynb):
            jupynb = self._dsource
            logger.info(
                f"Added notebook {jupynb.name!r} ({jupynb.id}, {jupynb.v}) by"
                f" user {settings.user.handle}."
            )

        self._print_logging_table()

        self._committed = True
        if isinstance(self._dsource, core.jupynb):
            finalize_publish(version=jupynb_v, calling_statement="commit(")

    def _print_logging_table(
        self, message: str = "Ingested the following dobjects:"
    ) -> None:
        """Pretty print logging messages."""
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
