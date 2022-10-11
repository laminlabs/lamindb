from pathlib import Path
from typing import Any, Dict, List, Union

import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from nbproject import dev, meta

from ..dev.db import Staged
from ..dev.db._query import query
from ..schema import core


class Ingest:
    """Ingest data.

    Ingest is an operation that stores and annotates data objects.

    Guide: :doc:`/db/guide/ingest`.

    Args:
        dsource: A data source. If `None` assumes a Jupyter Notebook. In the
            default schema `jupynb` and `pipeline_run` are the two allowed data
            sources.
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
