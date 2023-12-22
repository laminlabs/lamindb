from typing import Type

from django.db.models import Max, Subquery
from lnschema_core import Artifact, Dataset, Registry
from lnschema_core.types import VisibilityChoice

from lamindb._query_set import QuerySet


def filter(Registry: Type[Registry], **expressions) -> QuerySet:
    """See :meth:`~lamindb.dev.Registry.filter`."""
    if Registry in {Artifact, Dataset}:
        # visibility is set to 0 unless expressions contains id or uid equality
        if not ("id" in expressions or "uid" in expressions):
            visibility = "visibility"
            if not any(e.startswith(visibility) for e in expressions):
                expressions[
                    visibility
                ] = VisibilityChoice.default.value  # default visibility
            # if visibility is None, do not apply a filter
            # otherwise, it would mean filtering for NULL values, which doesn't make
            # sense for a non-NULLABLE column
            elif visibility in expressions and expressions[visibility] is None:
                expressions.pop(visibility)
    qs = QuerySet(model=Registry)
    if len(expressions) > 0:
        if expressions.get("version") == "latest" and hasattr(
            qs.model, "initial_version"
        ):
            versionless_expressions = expressions.copy()
            versionless_expressions.pop("version")
            qs = qs.filter(**versionless_expressions)
            qs = _filter_query_set_by_latest_version(qs)
            return qs
        else:
            return qs.filter(**expressions)
    else:
        return qs


def filter_version_family(record: Registry) -> QuerySet:
    qs = QuerySet(model=record._meta.model)

    # Get records with the same initial_version_id
    qs = qs.filter(initial_version_id=record.initial_version_id)

    return qs


def _filter_query_set_by_latest_version(query_set):
    # First, we create a subquery that gets the latest created_at for each initial_version_id
    latest_dates = (
        query_set.values("initial_version_id")
        .annotate(latest_date=Max("created_at"))
        .values("latest_date")
    )

    # Then, we filter the main queryset to only include records that have a created_at that matches the latest_date
    query_set = query_set.filter(created_at__in=Subquery(latest_dates))

    return query_set
