from __future__ import annotations

from typing import TYPE_CHECKING

from lnschema_core import Artifact, Collection

from ._query_set import QuerySet, process_expressions

if TYPE_CHECKING:
    from lnschema_core import Record


def filter(registry: type[Record], **expressions) -> QuerySet:
    """See :meth:`~lamindb.core.Record.filter`."""
    _using_key = None
    if "_using_key" in expressions:
        _using_key = expressions.pop("_using_key")
    expressions = process_expressions(registry, expressions)
    qs = QuerySet(model=registry, using=_using_key)
    if len(expressions) > 0:
        return qs.filter(**expressions)
    else:
        return qs
