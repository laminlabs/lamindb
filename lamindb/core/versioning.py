from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from lamin_utils import logger
from lamin_utils._base62 import CHARSET_DEFAULT as BASE62_CHARS
from lamindb_setup.core.upath import LocalPathClasses, UPath
from lnschema_core import ids

if TYPE_CHECKING:
    from lnschema_core.models import IsVersioned


def message_update_key_in_version_family(
    *,
    suid: str,
    existing_key: str,
    registry: str,
    new_key: str,
) -> str:
    return f'Or update key "{existing_key}" in your existing family:\n\nln.{registry}.filter(uid__startswith="{suid}").update(key="{new_key}")'


def increment_base62(s: str) -> str:
    # we don't need to throw an error for zzzz because uids are enforced to be unique
    # on the db level and have an enforced maximum length
    value = sum(BASE62_CHARS.index(c) * (62**i) for i, c in enumerate(reversed(s)))
    value += 1
    result = ""
    while value:
        value, remainder = divmod(value, 62)
        result = BASE62_CHARS[remainder] + result
    return result.zfill(len(s))


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
    if version is None and previous_version is not None:
        version = bump_version(previous_version, bump_type="major")
    return version


def create_uid(
    *,
    version: str | None = None,
    n_full_id: int = 20,
    revises: IsVersioned | None = None,
) -> tuple[str, IsVersioned | None]:
    if revises is not None:
        if not revises.is_latest:
            # need one more request
            revises = revises.__class__.objects.get(
                is_latest=True, uid__startswith=revises.stem_uid
            )
            logger.warning(
                f"didn't pass the latest version in `revises`, retrieved it: {revises}"
            )
        suid = revises.stem_uid
        vuid = increment_base62(revises.uid[-4:])
    else:
        suid = ids.base62(n_full_id - 4)
        vuid = "0000"
    if version is not None:
        if not isinstance(version, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0.1', '1', '2',"
                " etc."
            )
        if revises is not None:
            if version == revises.version:
                raise ValueError(
                    f"Please increment the previous version: '{revises.version}'"
                )
    return suid + vuid, revises


def get_new_path_from_uid(old_path: UPath, old_uid: str, new_uid: str):
    if isinstance(old_path, LocalPathClasses):
        # for local path, the rename target must be full path
        new_path = old_path.as_posix().replace(old_uid, new_uid)
    else:
        # for cloud path, the rename target must be the last part of the path
        new_path = old_path.name.replace(old_uid, new_uid)
    return new_path


def process_revises(
    revises: IsVersioned | None,
    version: str | None,
    name: str | None,
    type: type[IsVersioned],
) -> tuple[str, str, str, IsVersioned | None]:
    if revises is not None and not isinstance(revises, type):
        raise TypeError(f"`revises` has to be of type `{type.__name__}`")
    uid, revises = create_uid(
        revises=revises, version=version, n_full_id=type._len_full_uid
    )
    if revises is not None:
        if name is None:
            name = revises.name
    return uid, version, name, revises
