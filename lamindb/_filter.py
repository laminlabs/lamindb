from typing import Type

from django.db.models import OuterRef, Subquery
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


def _filter_query_set_by_latest_version(query_set):
    # Get records with the most recent created_at for each initial_version id
    subquery = query_set.filter(initial_version=OuterRef("initial_version")).order_by(
        "-created_at"
    )
    queryset1 = query_set.filter(
        initial_version__isnull=False,
        created_at=Subquery(subquery.values("created_at")[:1]),
    )

    # Get records where initial_version is null and exclude those records whose
    # id is an initial_version in queryset2
    queryset2 = query_set.filter(initial_version__isnull=True).exclude(
        id__in=queryset1.values_list("initial_version", flat=True)
    )

    # Combine the querysets
    final_queryset = queryset2.union(queryset1)

    return final_queryset
