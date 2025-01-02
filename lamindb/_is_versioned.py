from __future__ import annotations

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup.core.upath import UPath

from lamindb.models import IsVersioned

from ._utils import attach_func_to_class_method
from .core.versioning import create_uid, get_new_path_from_uid


# docstring handled through attach_func_to_class_method
def _add_to_version_family(self, revises: IsVersioned, version: str | None = None):
    old_uid = self.uid
    new_uid, revises = create_uid(revises=revises, version=version)
    if self.__class__.__name__ == "Artifact" and self._key_is_virtual:
        old_path = self.path
        new_path = get_new_path_from_uid(
            old_path=old_path, old_uid=old_uid, new_uid=new_uid
        )
        new_path = UPath(old_path).rename(new_path)
        logger.success(f"updated path from {old_path} to {new_path}!")
    self.uid = new_uid
    self.version = version
    self.save()
    logger.success(f"updated uid from {old_uid} to {new_uid}!")


METHOD_NAMES = [
    "_add_to_version_family",
]

if ln_setup._TESTING:  # type: ignore
    from inspect import signature

    SIGS = {name: signature(getattr(IsVersioned, name)) for name in METHOD_NAMES}

for name in METHOD_NAMES:
    attach_func_to_class_method(name, IsVersioned, globals())
