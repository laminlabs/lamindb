from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from lamin_utils import logger
from lamindb_setup.core.upath import LocalPathClasses, UPath
from lnschema_core import ids

if TYPE_CHECKING:
    from lnschema_core.models import IsVersioned


def bump_version(
    version: str,
    bump_type: str = "minor",
    behavior: Literal["prompt", "error", "ignore"] = "error",
) -> str:
    """Bumps the version number by major or minor depending on the bump_type flag.

    Parameters:
    version (str): The current version in "MAJOR" or "MAJOR.MINOR" format.
    bump_type (str): The type of version bump, either 'major' or 'minor'.

    Returns:
    str: The new version string.
    """
    try:
        # Split the version into major and minor parts if possible
        parts = version.split(".")
        major = int(parts[0])
        minor = int(parts[1]) if len(parts) > 1 else 0

        if bump_type == "major":
            # Bump the major version and reset the minor version
            new_version = f"{major + 1}"
        elif bump_type == "minor":
            # Bump the minor version
            new_version = f"{major}.{minor + 1}"
        else:
            raise ValueError("bump_type must be 'major' or 'minor'")

    except (ValueError, IndexError):
        if behavior == "prompt":
            new_version = input(
                f"The current version is '{version}' - please type the new version: "
            )
        elif behavior == "error":
            raise ValueError(
                "Cannot auto-increment non-integer castable version, please provide"
                " manually"
            ) from None
        else:
            logger.warning("could not auto-increment version, fix '?' manually")
            new_version = "?"
    return new_version


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
        version = bump_version(previous_version, bump_type="major")
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
    type: type[IsVersioned],
) -> tuple[str, str, str]:
    if is_new_version_of is not None and not isinstance(is_new_version_of, type):
        raise TypeError(f"is_new_version_of has to be of type {type}")
    if is_new_version_of is None:
        uid = init_uid(version=version, n_full_id=type._len_full_uid)
    else:
        uid, version = get_uid_from_old_version(is_new_version_of, version)
        if name is None:
            name = is_new_version_of.name
    return uid, version, name
