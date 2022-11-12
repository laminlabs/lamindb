from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from lnschema_core import Jupynb, Run
from nbproject import meta

from .dev.db import Staged
from .dev.db._select import select


class Ingest:
    """Ingest data objects for storage (files, arrays, `AnnData`, `DataFrame`, etc.).

    Store and link data objects through :meth:`lamindb.Ingest.add` and
    :meth:`lamindb.Ingest.commit`.

    Args:
        run: A data source. If `None` assumes a Jupyter Notebook. In the
            current core schema, `Jupynb` and `Run` are the two allowed
            data sources.

    For each staged data object, `Ingest` takes care of:

    1. Adding a record of :class:`~lamindb.schema.DObject` and linking it against
       its data source (:class:`~lamindb.schema.Run`).
    2. Linking features (:class:`lamindb.dev.db.Staged.link`) or other metadata
       (:class:`lamindb.dev.db.Staged.link_features`).
    3. Storing the corresponding data object in the storage location
       (:class:`~lamindb.schema.Storage`, `settings.instance.storage`).

    Guide: :doc:`/db/guide/ingest`.

    Examples:

    >>> import lamindb as ln
    >>> import lamindb as lns
    >>> ingest = ln.Ingest()  # run in Jupyter notebook or provide run!
    >>> filepath = ln.dev.datasets.file_jpg_paradisi05()
    >>> staged = ingest.add(filepath)
    >>> ingest.commit()
    """

    def _init_run(self, dsource: Union[Jupynb, Run]):
        if isinstance(dsource, Run):
            run = select(Run).where(Run.id == dsource.id).one_or_none()
            if run is None:
                run = dsource
            log = dict(run=f"{dsource.name!r} ({dsource.id})")
        elif isinstance(dsource, Jupynb):
            from lamindb._nb import _run as run  # type: ignore

            log = dict(jupynb=f"{dsource.name!r} ({dsource.id}, {dsource.v})")
        return run, log

    def __init__(self, run: Union[Jupynb, Run, None] = None):
        dsource = run  # rename
        if dsource is None:
            from lamindb._nb import _jupynb as jupynb

            if jupynb is not None:
                dsource = jupynb
            else:
                raise RuntimeError("Please provide a data source.")
        self._dsource = dsource  # data source (run or jupynb)
        self._run, self._runlog = self._init_run(dsource)
        self._staged: Dict = {}  # staged dobjects
        self._logs: List = []  # logging messages
        self._userlog = dict(user=f"{settings.user.handle} ({settings.user.id})")
        self._committed = False

    @property
    def run(self):
        """Run that is the data source for ingested objects."""
        return self._run

    def add(
        self,
        data: Any,
        *,
        name: Optional[str] = None,
        dobject_id: Optional[str] = None,
        adata_format: Optional[str] = None,
    ) -> Staged:
        """Stage data object for ingestion.

        Returns a staged object that can be linked against metadata, see
        :class:`~lamindb.dev.db.Staged`.

        Args:
            data: filepath or in-memory objects
            name: name of the data object, required of an in-memory object is passed
            dobject_id: id of the dobject
            adata_format: Use `h5ad` or `zarr` to store an `AnnData` object
        """
        staged = Staged(
            data,
            run=self.run,
            name=name,
            dobject_id=dobject_id,
            adata_format=adata_format,
        )
        self._staged[staged._filepath.as_posix()] = staged
        return staged

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

    def commit(
        self,
        *,
        nb_v: str = None,
        i_confirm_i_saved: bool = False,
        use_fsspec: bool = True,
    ) -> None:
        """Commit data to object storage and database.

        Args:
            nb_v: Notebook version to publish.
            i_confirm_i_saved: Only relevant outside Jupyter Lab & Notebook as a
                safeguard against losing the editor buffer content.
            use_fsspec: use fsspec to upload files.
        """
        if self._committed:
            logger.error("Already committed")
            return None

        if isinstance(self._dsource, Jupynb):
            from nbproject._publish import (
                finalize_publish,
                run_checks_for_publish,
                set_version,
            )

            result = run_checks_for_publish(
                calling_statement="commit(", i_confirm_i_saved=i_confirm_i_saved
            )
            if result != "checks-passed":
                return result

            # version to be set in finalize_publish()
            self._dsource.v = set_version(nb_v)

            # in case the nb exists already, update that entry
            result = (
                select(Jupynb)
                .where(Jupynb.id == self._dsource.id, Jupynb.v == self._dsource.v)
                .one_or_none()
            )
            if result is not None:
                self._dsource = result
                self._dsource.name = meta.live.title

            # also update run
            self.run.jupynb_v = self._dsource.v

        # insert dsource and run
        with sqm.Session(settings.instance.db_engine()) as session:
            session.add(self._dsource)
            session.commit()  # to satisfy foreign key constraint
            session.add(self.run)
            session.commit()
            # need to refresh here so that the both objects
            # are available for downstream use
            session.refresh(self._dsource)
            session.refresh(self.run)
        # sync db after changing locally
        settings.instance._update_cloud_sqlite_file()

        # one run that commits all dobjects
        for filepath_str, staged in self._staged.items():
            # TODO: run the appropriate clean-up operations if any aspect
            # of the ingestion fails
            staged._commit_dobject(use_fsspec=use_fsspec)
            self._logs.append({**staged._datalog, **self._runlog, **self._userlog})

        # one run that commits all linked entries
        for filepath_str, staged in self._staged.items():
            staged._commit_entries()

        if isinstance(self._dsource, Jupynb):
            jupynb = self._dsource
            logger.info(
                f"Added notebook {jupynb.name!r} ({jupynb.id}, {jupynb.v}) by"
                f" user {settings.user.handle}."
            )

        self._print_logging_table()

        self._committed = True
        if isinstance(self._dsource, Jupynb):
            finalize_publish(version=self._dsource.v, calling_statement="commit(")

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
