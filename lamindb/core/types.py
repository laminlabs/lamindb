from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

from anndata import AnnData
from lamindb_setup.core.types import UPathStr

from lamindb.base.types import (
    FeatureDtype,
    FieldAttr,
    ListLike,
    StrField,
    TransformType,
)

MuData = TypeVar("MuData")
SpatialData = TypeVar("SpatialData")

ScverseDataStructures = AnnData | MuData | SpatialData
