from typing import Type

from lnschema_core import Registry

from lamindb._queryset import QuerySet


def filter(Registry: Type[Registry], **expressions) -> QuerySet:
    """See :meth:`~lamindb.dev.Registry.filter`."""
    qs = QuerySet(model=Registry)
    if len(expressions) > 0:
        return qs.filter(**expressions)
    else:
        return qs
