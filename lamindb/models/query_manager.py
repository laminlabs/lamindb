from __future__ import annotations

import re
from functools import reduce
from typing import TYPE_CHECKING, NamedTuple

from django.db.models import (
    IntegerField,
    Manager,
    Q,
    QuerySet,
    TextField,
    Value,
)
from django.db.models.functions import Cast, Coalesce
from django.db.models.lookups import (
    Contains,
    Exact,
    IContains,
    IExact,
    IRegex,
    IStartsWith,
    Regex,
    StartsWith,
)
from lamin_utils import logger
from lamin_utils._lookup import Lookup
from lamindb_setup.core._docs import doc_args

if TYPE_CHECKING:
    from ..base.types import StrField


def _search(
    cls,
    string: str,
    *,
    field: StrField | list[StrField] | None = None,
    limit: int | None = 20,
    case_sensitive: bool = False,
    truncate_string: bool = False,
) -> QuerySet:
    """Search.

    Args:
        string: The input string to match against the field ontology values.
        field: The field or fields to search. Search all string fields by default.
        limit: Maximum amount of top results to return.
        case_sensitive: Whether the match is case sensitive.

    Returns:
        A sorted `DataFrame` of search results with a score in column `score`.
        If `return_queryset` is `True`.  `QuerySet`.

    See Also:
        :meth:`~lamindb.models.SQLRecord.filter`
        :meth:`~lamindb.models.SQLRecord.lookup`

    Examples:
        >>> ulabels = ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name")
        >>> ln.save(ulabels)
        >>> ln.ULabel.search("ULabel2")
    """
    if string is None:
        raise ValueError("Cannot search for None value! Please pass a valid string.")

    input_queryset = (
        cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    )
    registry = input_queryset.model
    name_field = getattr(registry, "_name_field", "name")
    if field is None:
        fields = [
            field.name
            for field in registry._meta.fields
            if field.get_internal_type() in {"CharField", "TextField"}
        ]
    else:
        if not isinstance(field, list):
            fields_input = [field]
        else:
            fields_input = field
        fields = []
        for field in fields_input:
            if not isinstance(field, str):
                try:
                    fields.append(field.field.name)
                except AttributeError as error:
                    raise TypeError(
                        "Please pass a SQLRecord string field, e.g., `CellType.name`!"
                    ) from error
            else:
                fields.append(field)

    if truncate_string:
        if (len_string := len(string)) > 5:
            n_80_pct = int(len_string * 0.8)
            string = string[:n_80_pct]

    string = string.strip()
    string_escape = re.escape(string)

    exact_lookup = Exact if case_sensitive else IExact
    regex_lookup = Regex if case_sensitive else IRegex
    contains_lookup = Contains if case_sensitive else IContains

    ranks = []
    contains_filters = []
    for field in fields:
        field_expr = Coalesce(
            Cast(field, output_field=TextField()),
            Value(""),
            output_field=TextField(),
        )
        # exact rank
        exact_expr = exact_lookup(field_expr, string)
        exact_rank = Cast(exact_expr, output_field=IntegerField()) * 200
        ranks.append(exact_rank)
        # exact synonym
        synonym_expr = regex_lookup(field_expr, rf"(?:^|.*\|){string_escape}(?:\|.*|$)")
        synonym_rank = Cast(synonym_expr, output_field=IntegerField()) * 200
        ranks.append(synonym_rank)
        # match as sub-phrase
        sub_expr = regex_lookup(
            field_expr, rf"(?:^|.*[ \|\.,;:]){string_escape}(?:[ \|\.,;:].*|$)"
        )
        sub_rank = Cast(sub_expr, output_field=IntegerField()) * 10
        ranks.append(sub_rank)
        # startswith and avoid matching string with " " on the right
        # mostly for truncated
        startswith_expr = regex_lookup(
            field_expr, rf"(?:^|.*\|){string_escape}[^ ]*(?:\|.*|$)"
        )
        startswith_rank = Cast(startswith_expr, output_field=IntegerField()) * 8
        ranks.append(startswith_rank)
        # match as sub-phrase from the left, mostly for truncated
        right_expr = regex_lookup(field_expr, rf"(?:^|.*[ \|]){string_escape}.*")
        right_rank = Cast(right_expr, output_field=IntegerField()) * 2
        ranks.append(right_rank)
        # match as sub-phrase from the right
        left_expr = regex_lookup(field_expr, rf".*{string_escape}(?:$|[ \|\.,;:].*)")
        left_rank = Cast(left_expr, output_field=IntegerField()) * 2
        ranks.append(left_rank)
        # simple contains filter
        contains_expr = contains_lookup(field_expr, string)
        contains_filter = Q(contains_expr)
        contains_filters.append(contains_filter)
        # also rank by contains
        contains_rank = Cast(contains_expr, output_field=IntegerField())
        ranks.append(contains_rank)
        # additional rule for truncated strings
        # weight matches from the beginning of the string higher
        # sometimes whole words get truncated and startswith_expr is not enough
        if truncate_string and field == name_field:
            startswith_lookup = StartsWith if case_sensitive else IStartsWith
            name_startswith_expr = startswith_lookup(field_expr, string)
            name_startswith_rank = (
                Cast(name_startswith_expr, output_field=IntegerField()) * 2
            )
            ranks.append(name_startswith_rank)

    ranked_queryset = (
        input_queryset.filter(reduce(lambda a, b: a | b, contains_filters))
        .alias(rank=sum(ranks))
        .order_by("-rank")
    )

    return ranked_queryset[:limit]


def _lookup(
    cls,
    field: StrField | None = None,
    return_field: StrField | None = None,
    using_key: str | None = None,
) -> NamedTuple:
    """Return an auto-complete object for a field.

    Args:
        field: The field to look up the values for. Defaults to first string field.
        return_field: The field to return. If `None`, returns the whole record.

    Returns:
        A `NamedTuple` of lookup information of the field values with a
        dictionary converter.

    See Also:
        :meth:`~lamindb.models.SQLRecord.search`

    Examples:
        >>> import bionty as bt
        >>> bt.settings.organism = "human"
        >>> bt.Gene.from_source(symbol="ADGB-DT").save()
        >>> lookup = bt.Gene.lookup()
        >>> lookup.adgb_dt
        >>> lookup_dict = lookup.dict()
        >>> lookup_dict['ADGB-DT']
        >>> lookup_by_ensembl_id = bt.Gene.lookup(field="ensembl_gene_id")
        >>> genes.ensg00000002745
        >>> lookup_return_symbols = bt.Gene.lookup(field="ensembl_gene_id", return_field="symbol")
    """
    from .sqlrecord import get_name_field

    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    field = get_name_field(registry=queryset.model, field=field)

    return Lookup(
        records=queryset,
        values=[i.get(field) for i in queryset.values()],
        tuple_name=cls.__class__.__name__,
        prefix="ln",
    ).lookup(
        return_field=(
            get_name_field(registry=queryset.model, field=return_field)
            if return_field is not None
            else None
        )
    )


# this is the default (._default_manager and ._base_manager) for lamindb models
class QueryManager(Manager):
    """Manage queries through fields.

    See Also:

        :class:`lamindb.models.QuerySet`
        `django Manager <https://docs.djangoproject.com/en/4.2/topics/db/managers/>`__

    Examples:

        >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name"))  # noqa
        >>> labels = ln.ULabel.filter(name__icontains = "label").all()
        >>> ln.ULabel(name="ULabel1").save()
        >>> label = ln.ULabel.get(name="ULabel1")
        >>> label.parents.set(labels)
        >>> manager = label.parents
        >>> manager.df()
    """

    def _track_run_input_manager(self):
        if hasattr(self, "source_field_name") and hasattr(self, "target_field_name"):
            if (
                self.source_field_name == "collection"
                and self.target_field_name == "artifact"
            ):
                from lamindb import settings
                from lamindb.core._context import context
                from lamindb.models.artifact import (
                    WARNING_RUN_TRANSFORM,
                    _track_run_input,
                )

                if (
                    context.run is None
                    and not settings.creation.artifact_silence_missing_run_warning
                ):
                    logger.warning(WARNING_RUN_TRANSFORM)
                _track_run_input(self.instance)

    def list(self, field: str | None = None):
        """Populate a list with the results.

        Examples:
            >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name"))
            >>> labels = ln.ULabel.filter(name__icontains="label").all()
            >>> ln.ULabel(name="ULabel1").save()
            >>> label = ln.ULabel.get(name="ULabel1")
            >>> label.parents.set(labels)
            >>> label.parents.list()
            >>> label.parents.list("name")
            ['ULabel1', 'ULabel2', 'ULabel3']
        """
        if field is None:
            return list(self.all())
        else:
            self._track_run_input_manager()
            return list(self.values_list(field, flat=True))

    def df(self, **kwargs):
        """Convert to DataFrame.

        For `**kwargs`, see :meth:`lamindb.models.QuerySet.df`.
        """
        return self.all().df(**kwargs)

    def all(self):
        """Return QuerySet of all.

        For `**kwargs`, see :meth:`lamindb.models.QuerySet.df`.
        """
        self._track_run_input_manager()
        return super().all()

    @doc_args(_search.__doc__)
    def search(self, string: str, **kwargs):
        """{}"""  # noqa: D415
        return _search(cls=self.all(), string=string, **kwargs)

    @doc_args(_lookup.__doc__)
    def lookup(self, field: StrField | None = None, **kwargs) -> NamedTuple:
        """{}"""  # noqa: D415
        return _lookup(cls=self.all(), field=field, **kwargs)

    def get_queryset(self):
        from .query_set import BasicQuerySet

        # QueryManager returns BasicQuerySet because it is problematic to redefine .filter and .get
        # for a query set used by the default manager
        return BasicQuerySet(model=self.model, using=self._db, hints=self._hints)
