from __future__ import annotations

from django.db.models import QuerySet
from lnschema_core import Artifact, Collection, Record
from lnschema_core.models import IsVersioned
from lnschema_core.types import VisibilityChoice

from lamindb import settings


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


def get(
    registry_or_queryset: type[Record] | QuerySet,
    idlike: int | str | None = None,
    **expressions,
) -> Record:
    if isinstance(registry_or_queryset, QuerySet):
        qs = registry_or_queryset
        registry = qs.model
    else:
        qs = QuerySet(model=registry_or_queryset)
        registry = registry_or_queryset
    if isinstance(idlike, int):
        return qs.get(id=idlike)
    elif isinstance(idlike, str):
        qs = qs.get(uid__startswith=idlike)
        if issubclass(registry, IsVersioned):
            if len(idlike) <= registry._len_stem_uid:
                return qs.latest_version().one()
            else:
                return qs.one()
        else:
            return qs.one()
    else:
        assert idlike is None  # noqa: S101
        # below behaves exactly like `.one()`
        return registry.get(**expressions)
