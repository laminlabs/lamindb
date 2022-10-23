import lnschema_core as core
from lndb_setup import settings


def track_usage(dobject_id, usage_type: core.type.usage):
    from ._add import add

    usage = add(
        core.usage(
            type=usage_type,
            user_id=settings.user.id,
            dobject_id=dobject_id,
        )
    )

    return usage.id
