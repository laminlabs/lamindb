from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

from lamindb_setup.types import UPathStr

from lamindb.base.types import (
    Dtype,
    FieldAttr,
    ListLike,
    StrField,
    TransformKind,
)

MuData = TypeVar("MuData")
SpatialData = TypeVar("SpatialData")

if TYPE_CHECKING:
    from anndata import AnnData

    ScverseDataStructures = AnnData | MuData | SpatialData
else:
    ScverseDataStructures = (
        Any  # AnnData | MuData | SpatialData; lazy to avoid importing anndata
    )
