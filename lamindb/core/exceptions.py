import warnings

warnings.warn(
    "`lamindb.core.exceptions` is deprecated, use `lamindb.errors` instead.",
    DeprecationWarning,
    stacklevel=2,
)

from ..errors import *  # noqa: F403 backward compat
