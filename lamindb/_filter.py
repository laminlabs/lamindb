from typing import Type

from lnschema_core import ORM

from lamindb._queryset import QuerySet


def filter(ORM: Type[ORM], **expressions) -> QuerySet:
    """See :meth:`~lamindb.dev.ORM.filter`."""
    qs = QuerySet(model=ORM)
    if len(expressions) > 0:
        return qs.filter(**expressions)
    else:
        return qs
