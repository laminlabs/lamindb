import warnings

warnings.warn(
    "`lamindb.base.ids` is deprecated, use `lamindb.base.uids` instead.",
    DeprecationWarning,
    stacklevel=2,
)

from .uids import *  # noqa: F403
