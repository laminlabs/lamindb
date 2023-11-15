from typing import Type

from lnschema_core import Registry

from lamindb._query_set import QuerySet


def filter(Registry: Type[Registry], **expressions) -> QuerySet:
    """See :meth:`~lamindb.dev.Registry.filter`."""
    qs = QuerySet(model=Registry)
    return qs
