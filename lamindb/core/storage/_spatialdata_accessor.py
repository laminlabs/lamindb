from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

from ._anndata_accessor import AnnDataAccessor

if TYPE_CHECKING:
    from zarr import Group

    from lamindb import Artifact


class _TablesAccessor:
    def __init__(self, tables: Group, artifact: Artifact | None = None):
        self._tables = tables

        self._artifact = artifact

    def __getitem__(self, key: str) -> AnnDataAccessor:
        return AnnDataAccessor(
            connection=None,
            storage=self._tables[key],
            filename=key,
            artifact=self._artifact,
        )

    def keys(self) -> list[str]:
        return list(self._tables.keys())

    def __repr__(self) -> str:
        """Description of the _TablesAccessor object."""
        descr = (
            f"Accessor for the SpatialData attribute tables\n  with keys: {self.keys()}"
        )
        return descr


class SpatialDataAccessor:
    """Cloud-backed SpatialData.

    For now only allows to access `tables`.
    """

    def __init__(self, storage: Group, name: str, artifact: Artifact | None = None):
        self.storage = storage
        self._name = name

        self._artifact = artifact

    @cached_property
    def tables(self) -> _TablesAccessor:
        """tables of the underlying SpatialData object."""
        return _TablesAccessor(self.storage["tables"], self._artifact)

    def __repr__(self):
        """Description of the SpatialDataAccessor object."""
        descr = (
            "SpatialDataAccessor object"
            f"\n  constructed for the SpatialData object {self._name}"
            f"\n    with tables: {self.tables.keys()}"
        )
        return descr
