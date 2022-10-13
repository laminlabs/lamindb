from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from nbproject import dev, meta

from ..dev.db import Staged
from ..dev.db._select import select
from ..schema import core


def set_nb_version(version):
    if version is not None:
        return version
    else:
        if meta.store.version == "draft":
            version = "1"
        else:
            version = meta.store.version
        return version


class Ingest:
    """Ingest data.

    Ingest is an operation that stores and annotates data objects.

    Guide: :doc:`/db/guide/ingest`.

    Args:
        dsource: A data source. If `None` assumes a Jupyter Notebook. In the
            default schema `jupynb` and `pipeline_run` are the two allowed data
            sources.
    """

    def _init_dtransform(self, dsource: Union[core.jupynb, core.pipeline_run]):
        if isinstance(dsource, core.pipeline_run):
            dtransform = select.dtransform(  # type: ignore
                pipeline_run_id=dsource.id
            ).one_or_none()
            if dtransform is None:
                dtransform = core.dtransform(pipeline_run_id=dsource.id)
            log = dict(pipeline_run=f"{dsource.name!r} ({dsource.id})")
        elif isinstance(dsource, core.jupynb):
            dtransform = select.dtransform(  # type: ignore
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
            dtransform=self._dtransform,
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

    def commit(self, *, nb_v: str = None, i_confirm_i_saved: bool = False) -> None:
        """Commit data to object storage and database.

        Args:
            nb_v: Notebook version to publish. Is automatically bumped from
                "draft" to "1" if `None`.
            i_confirm_i_saved: Only relevant outside Jupyter Lab & Notebook as a
                safeguard against losing the editor buffer content.
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
            self._dsource.v = set_nb_version(nb_v)

            # in case the nb exists already, update that entry
            result = select.jupynb(id=self._dsource.id, v=self._dsource.v).one_or_none()  # type: ignore  # noqa
            if result is not None:
                self._dsource = result
                self._dsource.name = meta.live.title

            # also update dtransform
            self._dtransform.jupynb_v = self._dsource.v

        # insert dsource and dtransform
        with sqm.Session(settings.instance.db_engine()) as session:
            session.add(self._dsource)
            session.commit()  # to satisfy foreign key constraint
            session.add(self._dtransform)
            session.commit()
            # need to refresh here so that the both objects
            # are available for downstream use
            session.refresh(self._dsource)
            session.refresh(self._dtransform)
        # sync db after changing locally
        settings.instance._update_cloud_sqlite_file()

        # one run that commits all dobjects
        for filepath_str, staged in self._staged.items():
            # TODO: run the appropriate clean-up operations if any aspect
            # of the ingestion fails
            staged._commit_dobject()
            self._logs.append(
                {**staged._datalog, **self._dtransformlog, **self._userlog}
            )

        # one run that commits all linked entries
        for filepath_str, staged in self._staged.items():
            staged._commit_entries()

        if isinstance(self._dsource, core.jupynb):
            jupynb = self._dsource
            logger.info(
                f"Added notebook {jupynb.name!r} ({jupynb.id}, {jupynb.v}) by"
                f" user {settings.user.handle}."
            )

        self._print_logging_table()

        self._committed = True
        if isinstance(self._dsource, core.jupynb):
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
