from __future__ import annotations

from lnschema_core import Artifact, Collection, Record
from lnschema_core.types import VisibilityChoice

from lamindb import settings

from ._query_set import QuerySet


def filter(registry: type[Record], **expressions) -> QuerySet:
    """See :meth:`~lamindb.core.Record.filter`."""
    _using_key = None
    if "_using_key" in expressions:
        _using_key = expressions.pop("_using_key")
    if registry in {Artifact, Collection}:
        # visibility is set to 0 unless expressions contains id or uid equality
        if not (
            "id" in expressions
            or "uid" in expressions
            or "uid__startswith" in expressions
        ):
            visibility = "visibility"
            if not any(e.startswith(visibility) for e in expressions):
                expressions[visibility] = (
                    VisibilityChoice.default.value
                )  # default visibility
            # if visibility is None, do not apply a filter
            # otherwise, it would mean filtering for NULL values, which doesn't make
            # sense for a non-NULLABLE column
            elif visibility in expressions and expressions[visibility] is None:
                expressions.pop(visibility)
    qs = QuerySet(model=registry, using=_using_key)
    if len(expressions) > 0:
        return qs.filter(**expressions)
    else:
        return qs
