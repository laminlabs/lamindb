"""Base types.

Classes
-------

.. autoclass:: CanonicalSuffix

Simple types
------------

.. autoclass:: ArtifactKind
.. autoclass:: TransformKind
.. autoclass:: BlockKind
.. autoclass:: BranchStatus
.. autoclass:: RunStatus
.. autoclass:: SimpleDtype
.. autoclass:: SimpleDtypeStr
.. autoclass:: SimpleDvalue
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
from collections.abc import Sequence
from typing import Literal, Union

from django.db.models.query_utils import DeferredAttribute as FieldAttr
from lamindb_setup.core.canonical_suffix import CanonicalSuffix  # noqa: F401
from lamindb_setup.types import AnyPathStr  # noqa: F401

# need to use Union because __future__.annotations doesn't do the job here <3.10
# typing.TypeAlias, >3.10 on but already deprecated
# Examples of list-like inputs: list[str], pandas.Series, numpy.ndarray.
ListLike = Sequence[str]
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

SimpleDvalue = int | float | str | bool | datetime.date | datetime.datetime | dict
"""Values corresponding to :class:`~lamindb.base.types.SimpleDtype`."""

SimpleDtype = (
    type[int]
    | type[float]
    | type[str]
    | type[bool]
    | type[datetime.date]
    | type[datetime.datetime]
    | type[dict]
)
"""Python types for simple scalar dtypes.

See section :ref:`Data types <dtypes-note>` on the :class:`~lamindb.Feature` page for more background.
"""

SimpleDtypeStr = Literal[
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
"""String-serialized representations for :class:`~lamindb.base.types.SimpleDtype`."""
DtypeStr = SimpleDtypeStr  # backward compat
Dtype = DtypeStr  # backward compat
DtypeObject = SimpleDvalue  # backward compat

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
