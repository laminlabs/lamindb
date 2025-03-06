from __future__ import annotations

from typing import TYPE_CHECKING

from lamindb_setup.core.types import UPathStr

from lamindb.base.types import (
    FeatureDtype,
    FieldAttr,
    ListLike,
    StrField,
    TransformType,
)
from lamindb.core._compat import is_package_installed

if TYPE_CHECKING:
    from anndata import AnnData
    from mudata import MuData
    from spatialdata import SpatialData

from anndata import AnnData

if is_package_installed("mudata"):
    from mudata import MuData
else:
    MuData = type("MuData", (), {})

if is_package_installed("spatialdata"):
    from spatialdata import SpatialData
else:
    SpatialData = type("SpatialData", (), {})

ScverseDataStructures = AnnData | MuData | SpatialData
