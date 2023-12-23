from typing import TYPE_CHECKING, Optional

from lnschema_core.models import TRANSFORM_TYPE_DEFAULT, Transform

from .dev.versioning import get_uid_from_old_version, init_uid

if TYPE_CHECKING:
    from lnschema_core.types import TransformType


def __init__(transform: Transform, *args, **kwargs):
    if len(args) == len(transform._meta.concrete_fields):
        super(Transform, transform).__init__(*args, **kwargs)
        return None
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    short_name: Optional[str] = (
        kwargs.pop("short_name") if "short_name" in kwargs else None
    )
    is_new_version_of: Optional[Transform] = (
        kwargs.pop("is_new_version_of") if "is_new_version_of" in kwargs else None
    )
    (kwargs.pop("initial_version_id") if "initial_version_id" in kwargs else None)
    version: Optional[str] = kwargs.pop("version") if "version" in kwargs else None
    type: Optional[TransformType] = (
        kwargs.pop("type") if "type" in kwargs else TRANSFORM_TYPE_DEFAULT
    )
    reference: Optional[str] = (
        kwargs.pop("reference") if "reference" in kwargs else None
    )
    # below is internal use that we'll hopefully be able to eliminate
    uid: Optional[str] = kwargs.pop("uid") if "uid" in kwargs else None
    if not len(kwargs) == 0:
        raise ValueError(
            "Only name, short_name, version, type, is_new_version_of can be passed,"
            f" but you passed: {kwargs}"
        )
    if is_new_version_of is None:
        new_uid = init_uid(version=version, n_full_id=Transform._len_full_uid)
    else:
        if not isinstance(is_new_version_of, Transform):
            raise TypeError("is_new_version_of has to be of type ln.Transform")
        new_uid, version = get_uid_from_old_version(
            is_new_version_of, version, n_full_id=Transform._len_full_uid
        )
        if name is None:
            name = is_new_version_of.name

    # this is only because the user-facing constructor allows passing an id
    # most others don't
    if uid is None:
        has_consciously_provided_uid = False
        uid = new_uid
    else:
        has_consciously_provided_uid = True
    super(Transform, transform).__init__(
        uid=uid,
        name=name,
        short_name=short_name,
        type=type,
        version=version,
        reference=reference,
        _has_consciously_provided_uid=has_consciously_provided_uid,
    )


Transform.__init__ = __init__
