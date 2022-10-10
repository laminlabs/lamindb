from pathlib import Path
from typing import Any, Optional, Union

from lndb_setup import settings
from nbproject import dev, meta

from ..._logger import logger
from ...dev import get_name_suffix_from_filepath, track_usage
from ...dev.file import store_file
from ...dev.object import infer_suffix, write_to_file
from ...schema import core
from .._insert import insert
from ._link_ingest import LinkIngest


class init_ingest:
    @classmethod
    def jupynb(cls) -> Optional[core.jupynb]:
        """Check if currently inside a Jupyter Notebook."""
        if dev.notebook_path() is not None:
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
        self._dobject = core.dobject(name=name, suffix=suffix)
        self._dobject.id = dobject_id if dobject_id is not None else self.dobject.id

        # access to the feature model
        self._feature_model = None  # feature model

        # access to the link operations
        self._link = LinkIngest(self)

        # dtransform
        self._dtransform = None

    @property
    def data(self) -> Any:
        """Data provided by the user upon init."""
        return self._data

    @property
    def dobject(self) -> core.dobject:
        """An dobject entry to be inserted."""
        return self._dobject

    @property
    def link(self) -> LinkIngest:
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

    def cleanup(self) -> None:
        """Clean up all linked entries."""
        if "dtransform" in self._link._entries:
            self._dtransform = None
        self._link._entries = {}

    def commit(self) -> None:
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
        insert.dobject_from_dtransform(  # type:ignore
            dobject=self.dobject, dtransform_id=self.dtransform.id  # type:ignore
        )

        # insert features and link to dobject
        if self.feature_model is not None:
            self.feature_model["model"].ingest(
                self.dobject.id, self.feature_model["df_curated"]
            )

        track_usage(self.dobject.id, usage_type="ingest")


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


def remove(filepath: Union[str, Path]) -> None:
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


def commit(jupynb_v: str = None, i_confirm_i_saved: bool = False) -> None:
    """Complete ingestion.

    Args:
        jupynb_v: Notebook version to publish. Is automatically set if `None`.
        i_confirm_i_saved: Only relevant outside Jupyter Lab as a safeguard against
            losing the editor buffer content because of accidentally publishing.
    """
    if jupynb is None:
        raise NotImplementedError
    else:
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


def print_logging_table(message: str = "Ingested the following dobjects:") -> None:
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
