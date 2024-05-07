from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Tuple

from lamindb_setup.core.upath import LocalPathClasses, UPath
from lnschema_core import ids

if TYPE_CHECKING:
    from lnschema_core.models import IsVersioned


def set_version(version: str | None = None, previous_version: str | None = None):
    """(Auto-) set version.

    If `version` is `None`, returns the stored version.
    Otherwise sets the version to the passed version.

    Args:
        version: Version string.
        previous_version: Previous version string.
    """
    if version == previous_version:
        raise ValueError(f"Please increment the previous version: '{previous_version}'")
    if version is None and previous_version is not None:
        try:
            version = str(int(previous_version) + 1)  # increment version by 1
        except ValueError:
            raise ValueError(
                "Cannot auto-increment non-integer castable version, please provide"
                " manually"
            ) from None
    return version


# uses `initial_version_id` to extract a stem_id that's part of id
def init_uid(
    *,
    version: str | None = None,
    n_full_id: int = 20,
    is_new_version_of: IsVersioned | None = None,
) -> str:
    if is_new_version_of is not None:
        stem_uid = is_new_version_of.stem_uid
    else:
        stem_uid = ids.base62(n_full_id - 4)
    if version is not None:
        if not isinstance(version, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0.1', '1', '2',"
                " etc."
            )
    return stem_uid + ids.base62_4()


def get_uid_from_old_version(
    is_new_version_of: IsVersioned,
    version: str | None = None,
    using_key: str | None = None,
) -> tuple[str, str]:
    """{}."""
    msg = ""
    if is_new_version_of.version is None:
        previous_version = "1"
        msg = f"setting previous version to '{previous_version}'"
    else:
        previous_version = is_new_version_of.version
    version = set_version(version, previous_version)
    new_uid = init_uid(
        version=version,
        n_full_id=is_new_version_of._len_full_uid,
        is_new_version_of=is_new_version_of,
    )
    # the following covers the edge case where the old file was unversioned
    if is_new_version_of.version is None:
        is_new_version_of.version = previous_version
        is_new_version_of.save(using=using_key)
        if msg != "":
            msg += f"& new version to '{version}'"
    return new_uid, version


def get_new_path_from_uid(old_path: UPath, old_uid: str, new_uid: str):
    if isinstance(old_path, LocalPathClasses):
        # for local path, the rename target must be full path
        new_path = old_path.as_posix().replace(old_uid, new_uid)
    else:
        # for cloud path, the rename target must be the last part of the path
        new_path = old_path.name.replace(old_uid, new_uid)
    return new_path


def process_is_new_version_of(
    is_new_version_of: IsVersioned,
    version: str | None,
    name: str | None,
    n_full_id: int,
) -> tuple[str, str, str]:
    if is_new_version_of is None:
        uid = init_uid(version=version, n_full_id=n_full_id)
    else:
        uid, version = get_uid_from_old_version(is_new_version_of, version)
        if name is None:
            name = is_new_version_of.name
    return uid, version, name
