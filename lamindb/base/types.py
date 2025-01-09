"""Types.

Central object types.

.. autosummary::
   :toctree: .

   ArtifactKind
   TransformType
   FeatureDtype

Basic types.

.. autosummary::
   :toctree: .

   UPathStr
   StrField
   ListLike
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

TransformType = Literal["pipeline", "notebook", "upload", "script", "function", "glue"]
ArtifactKind = Literal["dataset", "model"]
FeatureDtype = Literal[
    "cat",  # categorical variables
    "num",  # numerical variables
    "str",  # string variables
    "int",  # integer variables
    "float",  # float variables
    "bool",  # boolean variables
    "date",  # date variables
    "datetime",  # datetime variables
    "object",  # this is a pandas type, we're only using it for complicated types, not for strings
]
