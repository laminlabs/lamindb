"""Types.

.. autosummary::
   :toctree: .

   PathLike
   DataLike
   TransformType
"""
from pathlib import Path
from typing import Any, TypeVar

from upath import UPath

PathLike = TypeVar("PathLike", str, Path, UPath)
# statically typing the following is hard because these are all heavy
# dependencies, even DataFrame is heavy & slow to import
DataLike = Any

from lnschema_core.types import TransformType  # noqa
