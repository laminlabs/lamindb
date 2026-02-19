"""Base types.

Central object types
--------------------

.. autoclass:: ArtifactKind
.. autoclass:: TransformKind
.. autoclass:: BlockKind
.. autoclass:: BranchStatus
.. autoclass:: RunStatus
.. autoclass:: DtypeStr

Basic types
-----------

.. autoclass:: UPathStr
.. autoclass:: StrField
.. autoclass:: ListLike
.. autoclass:: FieldAttr
"""

from __future__ import annotations

import datetime
from typing import TYPE_CHECKING, Literal, Union

import numpy as np
from django.db.models.query_utils import DeferredAttribute as FieldAttr
from lamindb_setup.types import UPathStr  # noqa: F401

if TYPE_CHECKING:
    import pandas as pd

# need to use Union because __future__.annotations doesn't do the job here <3.10
# typing.TypeAlias, >3.10 on but already deprecated
# pd.Series as string to avoid importing pandas at runtime
ListLike = Union[list[str], "pd.Series", np.ndarray]
StrField = Union[str, FieldAttr]  # typing.TypeAlias

TransformKind = Literal["pipeline", "notebook", "script", "function"]
TransformType = TransformKind  # backward compat
ArtifactKind = Literal[
    "dataset", "model", "plan", "__lamindb_run__", "__lamindb_config__"
]
BlockKind = Literal["readme", "comment"]
"""Block kind, a `README.md`-type page or comment.

Any block expects Markdown as the formatting language.
"""

BranchStatus = Literal["standalone", "draft", "review", "merged", "closed"]
"""Branch status.

- `standalone`: Branch has no Merge Request intent.
- `draft`: Merge Request exists but is not ready for review.
- `review`: Merge Request is ready for review.
- `merged`: Merge Request has been merged into another branch.
- `closed`: Merge Request was closed without merging.
"""

RunStatus = Literal[
    "scheduled", "restarted", "started", "completed", "errored", "aborted"
]
"""Run status.

- `scheduled`: Run is scheduled (not yet started).
- `restarted`: Run was restarted.
- `started`: Run has started.
- `completed`: Run completed successfully.
- `errored`: Run ended with an error.
- `aborted`: Run was aborted.
"""

DtypeObject = int | float | str | bool | datetime.date | datetime.datetime | dict

DtypeStr = Literal[
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
"""String-serialized data type.

String-serialized representations of common data types.

Overview
========

============  ============  =================================================
description   lamindb       pandas
============  ============  =================================================
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

`lamindb` allows you to define a registry to which categoricals values are restricted.

For example, `'cat[ULabel]'` or `'cat[bionty.CellType]'` indicate that permissible values are stored in the `name` field of the `ULabel` or `CellType` registry, respectively.

You can also restrict to sub-types defined in registries via the `type` column, e.g., `'cat[ULabel[123456ABCDEFG]]'` indicates that values must be of the type with `uid="123456ABCDEFG"` within the `ULabel` registry.

Literal
=======

A `DtypeStr` object in `lamindb` is a `Literal` up to further specification of `"cat"`.

"""
Dtype = DtypeStr  # backward compat

RegistryId = Literal[
    "__lamindb_artifact__",
    "__lamindb_block__",
    "__lamindb_collection__",
    "__lamindb_feature__",
    "__lamindb_jsonvalue__",
    "__lamindb_project__",
    "__lamindb_record__",
    "__lamindb_run__",
    "__lamindb_schema__",
    "__lamindb_storage__",
    "__lamindb_transform__",
    "__lamindb_ulabel__",
]
