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
def init_id(
    *,
    provisional_id: Optional[None] = None,
    initial_version_id: Optional[str] = None,
    version: Optional[str] = None,
    n_full_id: int = 20,
) -> str:
    n_stem_id = n_full_id - 2
    if n_full_id == 20:
        gen_full_id = ids.base62_20
        gen_stem_id = ids.base62_18
    elif n_full_id:
        gen_full_id = ids.base62_14
        gen_stem_id = ids.base62_12
    if version is not None:
        if not isinstance(version, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0.1', '1', '2',"
                " etc."
            )
        # considered below, but not doing this right now
        # if version == "0":
        #     raise ValueError(
        #         "Please choose a version != '0', as it could be interpreted as `None`"
        #     )
    if initial_version_id is not None:
        stem_id = initial_version_id[:n_stem_id]
    else:
        stem_id = None
    # first consider an unversioned record
    if version is None and stem_id is None:
        provisional_id = gen_full_id()
        return provisional_id  # type: ignore
    # now consider a versioned record
    id_ext = ids.base62(2)
    if provisional_id is None and stem_id is None:
        stem_id = gen_stem_id()
        provisional_id = stem_id + id_ext
    elif stem_id is not None:
        assert isinstance(stem_id, str) and len(stem_id) == n_stem_id
        provisional_id = stem_id + id_ext
    return provisional_id  # type: ignore


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
) -> Tuple[str, str, str]:
    """{}"""
    msg = ""
    if is_new_version_of.version is None:
        previous_version = "1"
        msg = f"setting previous version to '{previous_version}'"
    else:
        previous_version = is_new_version_of.version
    version = set_version(version, previous_version)
    initial_version_id = get_initial_version_id(is_new_version_of)
    new_file_id = init_id(
        provisional_id=is_new_version_of.id,
        initial_version_id=initial_version_id,
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
    return new_file_id, initial_version_id, version  # type: ignore
