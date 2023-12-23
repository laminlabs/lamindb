from typing import Optional, Tuple, Union

from lnschema_core import ids
from lnschema_core.models import IsVersioned


def set_version(version: Optional[str] = None, previous_version: Optional[str] = None):
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
    version: Optional[str] = None,
    n_full_id: int = 20,
    is_new_version_of: Optional[IsVersioned] = None,
) -> str:
    if is_new_version_of is not None:
        stem_uid = is_new_version_of.stem_uid
    else:
        if n_full_id == 20:
            stem_uid = ids.base62_16()
        elif n_full_id == 16:
            stem_uid = ids.base62_12()
    if version is not None:
        if not isinstance(version, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0.1', '1', '2',"
                " etc."
            )
    return stem_uid + ids.base62_4()


def get_uid_from_old_version(
    is_new_version_of: IsVersioned,
    version: Optional[str],
    n_full_id: int = 20,
) -> Tuple[str, str]:
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
        n_full_id=n_full_id,
        is_new_version_of=is_new_version_of,
    )
    # the following covers the edge case where the old file was unversioned
    if is_new_version_of.version is None:
        is_new_version_of.version = previous_version
        is_new_version_of.save()
        if msg != "":
            msg += f"& new version to '{version}'"
    return new_uid, version
