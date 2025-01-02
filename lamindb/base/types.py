from __future__ import annotations

from enum import Enum
from typing import Literal, Union

import numpy as np
import pandas as pd
from django.db.models import IntegerChoices  # needed elsewhere
from django.db.models.query_utils import DeferredAttribute as FieldAttr

# need to use Union because __future__.annotations doesn't do the job here <3.10
# typing.TypeAlias, >3.10 on but already deprecated
ListLike = Union[list[str], pd.Series, np.array]
StrField = Union[str, FieldAttr]  # typing.TypeAlias

TransformType = Literal["pipeline", "notebook", "upload", "script", "function", "glue"]
ArtifactType = Literal["dataset", "model"]
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


class VisibilityChoice(IntegerChoices):
    default = 1
    hidden = 0
    trash = -1
