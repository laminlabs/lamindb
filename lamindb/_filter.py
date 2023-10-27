from typing import Type

from lnschema_core import Registry

from lamindb._query_set import QuerySet


def filter(Registry: Type[Registry], **expressions) -> QuerySet:
    """See :meth:`~lamindb.dev.Registry.filter`."""
    if Registry.__name__ in {"File", "Dataset"}:
        # visibility is set to <2 by default
        if not any([e.startswith("visibility") for e in expressions]):
            expressions["visibility__lt"] = 1
        # if visibility is None, will not apply any filter for visibility
        elif "visibility" in expressions:
            if expressions["visibility"] is None:
                expressions.pop("visibility")
            elif expressions["visibility"] == "default":
                expressions.pop("visibility")
                expressions["visibility__lt"] = 1
    qs = QuerySet(model=Registry, using=using)
    if len(expressions) > 0:
        return qs.filter(**expressions)
    else:
        return qs
