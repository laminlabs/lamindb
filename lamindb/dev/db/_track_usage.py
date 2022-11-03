from lndb_setup import settings
from lnschema_core import Usage
from lnschema_core.dev.type import usage as UsageType


def track_usage(dobject_id, usage_type: UsageType):
    from ._add import add

    usage = Usage(
        type=usage_type,
        user_id=settings.user.id,
        dobject_id=dobject_id,
    )
    add(usage)

    return usage.id
