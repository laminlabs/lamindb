import hashlib
from pathlib import Path
from sys import getsizeof
from typing import Any, Optional

import lnschema_core as core
from lndb_setup import settings

from .._core import get_name_suffix_from_filepath
from ..file import store_file, write_adata_zarr
from ..object import infer_suffix, write_to_file
from ._insert import insert
from ._linkstaged import LinkStaged
from ._track_usage import track_usage


class Staged:
    """Staged data objects, initiated upon :meth:`~lamindb.db.Ingest.add`.

    Guide: :doc:`/db/guide/ingest`.

    Args:
        data: filepath or in-memory objects
        dtransform: The data transform that links to the data source of the data object.
        name: name of the data object, required of an in-memory object is passed
        dobject_id: id of the dobject
        adata_format: Use `h5ad` or `zarr` to store an `AnnData` object
    """

    def __init__(
        self,
        data: Any,
        *,
        dtransform: core.dtransform,
        name: Optional[str] = None,
        dobject_id: Optional[str] = None,
        adata_format: Optional[str] = None,
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
            # if zarr, then filename is just name.zarr
            self._dmem = data
            suffix = infer_suffix(data, adata_format)
            if name is None:
                raise RuntimeError("Provide name if ingesting in-memory data.")
            self._filepath = Path(f"{name}{suffix}")
            # write to file
            # TODO: should this raise an error/warning if filepath exists?
            if suffix != ".zarr" and not self.filepath.exists():
                write_to_file(self.dmem, self.filepath)  # type: ignore

        self._dobject = core.dobject(name=name, suffix=suffix)
        self._dobject.id = dobject_id if dobject_id is not None else self.dobject.id
        # streamed
        if suffix != ".zarr":
            self._dobject.checksum = compute_checksum(self._filepath)

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

        See :class:`~lamindb.dev.db.LinkFeatureModel`
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
        dobject_storage_key = f"{self.dobject.id}{self.dobject.suffix}"

        if self.dobject.suffix != ".zarr":
            size = store_file(self.filepath, dobject_storage_key)
        else:
            # adata size
            size = getsizeof(self.dmem)
            storepath = settings.instance.storage.key_to_filepath(dobject_storage_key)
            write_adata_zarr(self.dmem, storepath)
        self._dobject.size = size

        # insert dobject first to satisfy foreign key constraints
        insert.dobject_from_dtransform(  # type:ignore
            dobject=self.dobject, dtransform_id=self._dtransform.id  # type:ignore
        )

        # insert all linked entries
        for table_name, entries in self.link.linked_entries.items():
            insert.from_list(table_name=table_name, entries=entries)  # type:ignore

        # insert features and link to dobject
        if self.feature_model is not None:
            self.feature_model["model"].ingest(
                self.dobject.id, self.feature_model["df_curated"]
            )

        track_usage(self.dobject.id, usage_type="ingest")


def compute_checksum(path: Path):
    # based on https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file  # noqa
    hash_md5 = hashlib.md5()
    with open(path, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
