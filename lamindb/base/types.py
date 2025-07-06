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
from lamindb_setup.types import UPathStr  # noqa: F401

# need to use Union because __future__.annotations doesn't do the job here <3.10
# typing.TypeAlias, >3.10 on but already deprecated
ListLike = Union[list[str], pd.Series, np.array]
StrField = Union[str, FieldAttr]  # typing.TypeAlias

TransformType = Literal[
    "pipeline", "notebook", "upload", "script", "function", "linker"
]
ArtifactKind = Literal["dataset", "model", "__lamindb_run__"]

# below is used for Feature.dtype and Param.dtype
Dtype = Literal[
    "cat",  # categoricals
    "num",  # numericals
    "str",  # string
    "int",  # integer / numpy.integer
    "float",  # float
    "bool",  # boolean
    "date",  # date
    "datetime",  # datetime
    "dict",  # dictionary
    "object",  # this is a pandas input dtype, we're only using it for complicated types, not for strings
    "path",  # path, validated as str, but specially treated in the UI
]
"""Data type.

String-serialized representations of common data types.

Overview
========

============  ============  =================================================
description   lamindb       pandas
============  ============  =================================================
categorical   `"cat"`       `category`
numerical     `"num"`       `int | float`
integer       `"int"`       `int64 | int32 | int16 | int8 | uint | ...`
float         `"float"`     `float64 | float32 | float16 | float8 | ...`
string        `"str"`       `object`
datetime      `"datetime"`  `datetime`
date          `"date"`      `object` (pandera requires an ISO-format string, convert with `df["date"] = df["date"].dt.date`)
dictionary    `"dict"`      `object`
path          `"path"`      `str` (pandas does not have a dedicated path type, validated as `str`)
============  ============  =================================================

Categoricals
============

Beyond indicating that a feature is a categorical, `lamindb` allows you to define the registry to which values are restricted.

For example, `'cat[ULabel]'` or `'cat[bionty.CellType]'` indicate that permissible values are from the `ULabel` or `CellType` registry, respectively.

You can also reference multiple registries, e.g., `'cat[ULabel|bionty.CellType]'` indicates that values can be from either registry.

You can also restrict to sub-types defined in registries via the `type` column, e.g., `'cat[ULabel[CellMedium]]'` indicates that values must be of type `CellMedium` within the `ULabel` registry.

Literal
=======

A `Dtype` object in `lamindb` is a `Literal` up to further specification of `"cat"`.

"""
FeatureDtype = Dtype  # backward compat
