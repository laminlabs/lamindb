"""Base types.

Central object types
--------------------

.. autoclass:: ArtifactKind
.. autoclass:: TransformKind
.. autoclass:: BlockKind
.. autoclass:: BranchStatus
.. autoclass:: RunStatus
.. autoclass:: SimpleDtype
.. autoclass:: DtypeStr

Basic types
-----------

.. autoclass:: AnyPathStr
.. autoclass:: StrField
.. autoclass:: ListLike
.. autoclass:: FieldAttr
"""

from __future__ import annotations

import datetime
from typing import TYPE_CHECKING, Literal, Union

import numpy as np
from django.db.models.query_utils import DeferredAttribute as FieldAttr
from lamindb_setup.types import AnyPathStr  # noqa: F401

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

=============  =====  ==================================================
status         code   description
=============  =====  ==================================================
`closed`       -2     Change Request was closed without merging.
`merged`       -1     The branch was merged into another branch.
`standalone`   0      A standalone branch without Change Request.
`draft`        1      Change Request exists but is not ready for review.
`review`       2      Change Request is ready for review.
=============  =====  ==================================================

The database stores the branch status as an integer code in field `_status_code`.
"""

RunStatus = Literal[
    "scheduled", "restarted", "started", "completed", "errored", "aborted"
]
"""Run status.

===========  =====  ===========================
status       code   description
===========  =====  ===========================
`scheduled`  -3     The run is scheduled.
`restarted`  -2     The run was restarted.
`started`    -1     The run has started.
`completed`  0      The run completed successfully.
`errored`    1      The run ended with an error.
`aborted`    2      The run was aborted.
===========  =====  ===========================

The database stores the run status as an integer code in field `_status_code`.
"""

RUN_STATUS_TO_CODE: dict[RunStatus, int] = {
    "scheduled": -3,
    "restarted": -2,
    "started": -1,
    "completed": 0,
    "errored": 1,
    "aborted": 2,
}
RUN_CODE_TO_STATUS: dict[int, RunStatus] = {
    code: status for status, code in RUN_STATUS_TO_CODE.items()
}

BRANCH_STATUS_TO_CODE: dict[BranchStatus, int] = {
    "closed": -2,
    "merged": -1,
    "standalone": 0,
    "draft": 1,
    "review": 2,
}
BRANCH_CODE_TO_STATUS: dict[int, BranchStatus] = {
    code: status for status, code in BRANCH_STATUS_TO_CODE.items()
}

DtypeObject = int | float | str | bool | datetime.date | datetime.datetime | dict

SimpleDtype = (
    type[int]
    | type[float]
    | type[str]
    | type[bool]
    | type[datetime.date]
    | type[datetime.datetime]
    | type[dict]
)
"""Native Python classes for LaminDB's simple scalar dtype categories.

This alias represents the preferred constructor inputs for simple feature dtypes
(`int`, `float`, `str`, `bool`, `datetime.date`, `datetime.datetime`, `dict`).
"""

DtypeStr = Literal[
    "num",  # numericals
    "int",  # integer / numpy.integer
    "float",  # float
    "str",  # string
    "bool",  # boolean
    "datetime",  # datetime
    "datetime64[ns, UTC]",  # timezone-aware datetime
    "date",  # date
    "dict",  # dictionary
    "path",  # path, validated as str, but specially treated in the UI
    "url",  # URL, validated as str, but specially treated in the UI
    "object",  # this is a pandas input dtype, we're only using it for complicated types, not for strings; consciously currently not documented
]
"""String-serialized representations of common data types.

===============  ====================  =================================================
description      lamindb (str)         pandas
===============  ====================  =================================================
numerical        `num`                 `int | float`
integer          `int`                 `int64 | int32 | int16 | int8 | uint | ...`
float            `float`               `float64 | float32 | float16 | float8 | ...`
string           `str`                 `object`
boolean          `bool`                `boolean | bool`
datetime (naive) `datetime`            `datetime`
datetime (tz)    `datetime64[ns, UTC]` `datetime64[ns, UTC]`
date             `date`                `object` (pandera requires an ISO-format string, convert with `df["date"] = df["date"].dt.date`)
dictionary       `dict`                `object`
path             `path`                `str` (pandas does not have a dedicated path type, validated as `str`)
url              `url`                 `str` (pandas does not have a dedicated url type, validated as `str`)
===============  ====================  =================================================

.. admonition:: Categorical and relational data types

    These are **not** contained in the `DTypeStr` `Literal`.

    For any categorical, you can restrict the permissible values to the values defined in a registry.
    When serializing this to a string, then `'cat[ULabel]'` or `'cat[bionty.CellType]'` indicate that permissible values are stored in the `name` field of the `ULabel` or `CellType` registry, respectively.
    You can also restrict to sub-types defined in registries via the `type` field, e.g., `'cat[ULabel[123456ABCDEFG]]'` indicates that values must be of the type with `uid="123456ABCDEFG"` within the `ULabel` registry.

    In LaminDB, categoricals define relationships with registries. See :class:`~lamindb.Feature` for more details.

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
