"""Types.

Central object types.

.. autosummary::
   :toctree: .

   ArtifactKind
   TransformType
   Dtype

Basic types.

.. autosummary::
   :toctree: .

   UPathStr
   StrField
   ListLike
   FieldAttr
"""

from __future__ import annotations

from typing import Literal, Union

import numpy as np
import pandas as pd
from django.db.models.query_utils import DeferredAttribute as FieldAttr
from lamindb_setup.core.types import UPathStr  # noqa: F401

# need to use Union because __future__.annotations doesn't do the job here <3.10
# typing.TypeAlias, >3.10 on but already deprecated
ListLike = Union[list[str], pd.Series, np.array]
StrField = Union[str, FieldAttr]  # typing.TypeAlias

TransformType = Literal[
    "pipeline", "notebook", "upload", "script", "function", "linker"
]
ArtifactKind = Literal["dataset", "model"]
Dtype = Literal[
    "cat",  # categoricals
    "num",  # numericals
    "str",  # string
    "int",  # integer / numpy.integer
    "float",  # float
    "bool",  # boolean
    "date",  # date
    "datetime",  # datetime
    "object",  # this is a pandas dtype, we're only using it for complicated types, not for strings
]
"""Data type.

============  =================================================
lamindb       numpy / pandas
============  =================================================
`"cat"`       `category`
`"num"`       `int | float`
`"int"`       `int64 | int32 | int16 | int8 | uint | ...`
`"float"`     `float64 | float32 | float16 | float8 | ...`
============  =================================================
"""
FeatureDtype = Dtype  # backward compat
