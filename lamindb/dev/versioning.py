from typing import Optional, Tuple, Union

from lnschema_core import ids
from lnschema_core.models import File, Transform


def set_version(version: Optional[str] = None, previous_version: Optional[str] = None):
    """(Auto-) set version.

    If `version` is `None`, returns the stored version.
    Otherwise sets the version to the passed version.

    Args:
        version: Version string.
        stored_version: Mock stored version for testing purposes.
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
            )
    return version


# uses `initial_version_id` to extract a stem_id that's part of id
def init_uid(
    *,
    version: Optional[str] = None,
    n_full_id: int = 20,
) -> str:
    if n_full_id == 20:
        gen_full_id = ids.base62_20
    elif n_full_id == 14:
        gen_full_id = ids.base62_14
    if version is not None:
        if not isinstance(version, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0.1', '1', '2',"
                " etc."
            )
    return gen_full_id()


def get_initial_version_id(is_new_version_of: Union[File, Transform]):
    if is_new_version_of.initial_version_id is None:
        initial_version_id = is_new_version_of.id
    else:
        initial_version_id = is_new_version_of.initial_version_id
    return initial_version_id


def get_ids_from_old_version(
    is_new_version_of: Union[File, Transform],
    version: Optional[str],
    n_full_id: int = 20,
) -> Tuple[str, int, str]:
    """{}"""
    msg = ""
    if is_new_version_of.version is None:
        previous_version = "1"
        msg = f"setting previous version to '{previous_version}'"
    else:
        previous_version = is_new_version_of.version
    version = set_version(version, previous_version)
    initial_version_id = get_initial_version_id(is_new_version_of)
    new_uid = init_uid(
        version=version,
        n_full_id=n_full_id,
    )
    # the following covers the edge case where the old file was unversioned
    if is_new_version_of.version is None:
        is_new_version_of.version = previous_version
        is_new_version_of.save()
        if msg != "":
            msg += (
                f"& new version to '{version}' (initial_version_id ="
                f" '{initial_version_id}')"
            )
    return new_uid, initial_version_id, version  # type: ignore
