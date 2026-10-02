from __future__ import annotations

import keyword
import re
from collections import namedtuple
from functools import reduce
from typing import TYPE_CHECKING, Any, Literal, NamedTuple

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
from lamindb_setup import logger
from lamindb_setup.core import deprecated
from lamindb_setup.core._docs import doc_args

if TYPE_CHECKING:
    from collections.abc import Iterable

    from ..base.types import StrField

SEARCH_QUERY_DEFAULT_LIMIT = 20


def _search(
    cls,
    string: str,
    *,
    field: StrField | list[StrField] | None = None,
    limit: int | None = SEARCH_QUERY_DEFAULT_LIMIT,
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

        ::

            records = ln.ULabel.from_values(["Label1", "Label2", "Label3"]).save()
            ln.ULabel.search("Label2")
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


def _append_records_to_list(df_dict: dict, value: str, record) -> None:
    """Append unique records to a list."""
    values_list = df_dict[value]

    if not isinstance(values_list, list):
        values_list = [values_list]
    try:
        df_dict[value] = list(dict.fromkeys(values_list + [record]))
    except TypeError:
        df_dict[value] = values_list


def _create_df_dict(
    df: Any = None,
    field: str | None = None,
    records: list | None = None,
    values: list | None = None,
    tuple_name: str | None = None,
) -> dict:
    """Create a dict with {lookup key: records in namedtuple}.

    Value is a list of namedtuples if multiple records match the same key.
    """
    if df is not None:
        records = df.itertuples(index=False, name=tuple_name)
        values = df[field]
    df_dict: dict = {}  # a dict of namedtuples as records and values as keys
    for i, row in enumerate(records):  # type:ignore
        value = values[i]  # type:ignore
        if not isinstance(value, str):
            continue
        if value == "":
            continue
        if value in df_dict:
            _append_records_to_list(df_dict=df_dict, value=value, record=row)
        else:
            df_dict[value] = row
    return df_dict


class _ListValueWrapper:
    """Wrapper that warns when a list value is accessed and applies keep strategy."""

    def __init__(
        self,
        field_name: str,
        values: list,
        keep: Literal["first", "last", False],
        return_field: str | None = None,
    ):
        self._field_name = field_name
        self._values = values
        self._keep = keep
        self._return_field = return_field
        self._accessed = False

    def _warn_and_process(self):
        """Issue warning and return processed value."""
        if not self._accessed:
            logger.warning(
                f"{len(self._values)} records found for '{self._field_name}'. "
                f"Returning based on keep='{self._keep}'."
            )
            self._accessed = True

        # Apply keep strategy
        if self._keep == "first":
            selected_value = self._values[0] if self._values else None
        elif self._keep == "last":
            selected_value = self._values[-1] if self._values else None
        elif self._keep is False:
            selected_value = self._values
        else:
            selected_value = self._values[0] if self._values else None

        # Apply return_field if specified
        if self._return_field is not None:
            if self._keep is False and isinstance(selected_value, list):
                return [
                    getattr(item, self._return_field)
                    if hasattr(item, self._return_field)
                    else item
                    for item in selected_value
                ]
            elif hasattr(selected_value, self._return_field):
                return getattr(selected_value, self._return_field)

        return selected_value

    def __getattr__(self, name):
        """Intercept any attribute access to trigger warning."""
        processed_value = self._warn_and_process()
        return getattr(processed_value, name)

    def __str__(self):
        """String representation triggers warning."""
        return str(self._warn_and_process())

    def __repr__(self):
        """Representation triggers warning."""
        return repr(self._warn_and_process())

    def __iter__(self):
        """Iteration triggers warning."""
        processed_value = self._warn_and_process()
        return iter(processed_value)

    def __len__(self):
        """Length check triggers warning."""
        processed_value = self._warn_and_process()
        return len(processed_value)

    def __getitem__(self, key):
        """Indexing triggers warning."""
        processed_value = self._warn_and_process()
        return processed_value[key]

    def __bool__(self):
        """Boolean conversion triggers warning."""
        processed_value = self._warn_and_process()
        return bool(processed_value)

    def __eq__(self, other):
        """Equality comparison triggers warning."""
        processed_value = self._warn_and_process()
        return processed_value == other


class Lookup:
    """Lookup object with dot and [] access."""

    # removed DataFrame type annotation to speed up import time
    def __init__(
        self,
        field: str | None = None,
        tuple_name="MyTuple",
        prefix: str = "bt",
        df: Any = None,
        values: Iterable | None = None,
        records: list | None = None,
        keep: Literal["first", "last", False] = "first",
    ) -> None:
        self._tuple_name = tuple_name
        if df is not None:
            if df.shape[0] > 500000:
                logger.warning(
                    "generating lookup object from >500k keys is not recommended and"
                    " extremely slow"
                )
            values = df[field]
        self._df_dict = _create_df_dict(
            df=df,
            field=field,
            records=records,
            values=values,  # type:ignore
            tuple_name=self._tuple_name,
        )
        lkeys = self._to_lookup_keys(values=values, prefix=prefix)  # type:ignore
        self._lookup_dict = self._create_lookup_dict(lkeys=lkeys, df_dict=self._df_dict)
        self._prefix = prefix
        self._keep = keep

    def _to_lookup_keys(self, values: Iterable, prefix: str) -> dict:
        """Convert a list of strings to tab-completion allowed formats.

        Returns:
            {lookup_key: value_or_values}
        """
        lkeys: dict = {}
        for value in list(values):
            if not isinstance(value, str):
                continue
            # replace any special character with _
            lkey = re.sub("[^0-9a-zA-Z_]+", "_", str(value)).lower()
            if lkey == "":  # empty strings are skipped
                continue
            if not lkey[0].isalpha():  # must start with a letter
                lkey = f"{prefix.lower()}_{lkey}"

            if lkey in lkeys:
                # if multiple values have the same lookup key
                # put the values into a list
                _append_records_to_list(df_dict=lkeys, value=lkey, record=value)
            else:
                lkeys[lkey] = value
        return lkeys

    def _create_lookup_dict(self, lkeys: dict, df_dict: dict) -> dict:
        lkey_dict: dict = {}  # a dict of namedtuples as records and lookup keys as keys
        for lkey, values in lkeys.items():
            if isinstance(values, list):
                combined_list = []
                for v in values:
                    records = df_dict.get(v)
                    if isinstance(records, list):
                        combined_list += records
                    else:
                        combined_list.append(records)
                lkey_dict[lkey] = combined_list
            else:
                lkey_dict[lkey] = df_dict.get(values)

        return lkey_dict

    def dict(self) -> dict:
        """Dictionary of the lookup."""
        return self._df_dict

    def lookup(self, return_field: str | None = None) -> NamedTuple:
        """Lookup records with dot access."""
        # Create a copy to avoid modifying the original
        lookup_dict_copy = self._lookup_dict.copy()

        # Process values, wrapping lists in warning wrapper
        processed_dict = {}
        duplicate_counts: dict[str, int] = {}
        for key, value in lookup_dict_copy.items():
            # Handle Python keywords by appending an underscore
            if keyword.iskeyword(key):
                key = f"{key}_"
            if isinstance(value, list) and len(value) > 1:
                if self._keep is False:
                    # Keep all duplicates and warn lazily on access.
                    processed_dict[key] = _ListValueWrapper(
                        key, value, self._keep, return_field
                    )
                else:
                    # Eagerly resolve duplicates for keep="first"/"last" so attribute
                    # values have the selected element type.
                    selected_value = value[0] if self._keep == "first" else value[-1]
                    duplicate_counts[key] = len(value)
                    if return_field is not None and hasattr(
                        selected_value, return_field
                    ):
                        processed_dict[key] = getattr(selected_value, return_field)
                    else:
                        processed_dict[key] = selected_value
            else:
                # Handle single values or single-item lists
                if isinstance(value, list) and len(value) == 1:
                    value = value[0]  # Unwrap single-item lists

                if return_field is not None and hasattr(value, return_field):
                    processed_dict[key] = getattr(value, return_field)
                else:
                    processed_dict[key] = value

        keys: list = list(processed_dict.keys()) + ["dict"]
        MyTuple = namedtuple("Lookup", keys)  # type:ignore

        if self._keep is not False and duplicate_counts:
            # Warn once per duplicated key on first attribute access while keeping
            # eager value resolution (non-wrapper return types).
            MyTuple._duplicate_counts = duplicate_counts  # type:ignore[attr-defined]
            MyTuple._warned_duplicate_keys = set()  # type:ignore[attr-defined]
            MyTuple._keep = self._keep  # type:ignore[attr-defined]

            def _lookup_getattribute(instance, name):
                cls = tuple.__getattribute__(instance, "__class__")
                duplicate_map = cls._duplicate_counts  # type:ignore[attr-defined]
                warned_keys = cls._warned_duplicate_keys  # type:ignore[attr-defined]
                if name in duplicate_map and name not in warned_keys:
                    logger.warning(
                        f"{duplicate_map[name]} records found for '{name}'. "
                        f"Returning based on keep='{cls._keep}'."  # type:ignore[attr-defined]
                    )
                    warned_keys.add(name)
                return tuple.__getattribute__(instance, name)

            MyTuple.__getattribute__ = _lookup_getattribute  # type:ignore[method-assign]

        return MyTuple(**processed_dict, dict=self.dict)  # type:ignore


def _lookup(
    cls,
    field: StrField | None = None,
    return_field: StrField | None = None,
    using: str | None = None,
    keep: Literal["first", "last", False] = "first",
) -> NamedTuple:
    """Return an auto-complete object for a field.

    Args:
        field: The field to look up the values for. Defaults to first string field.
        return_field: The field to return. If `None`, returns the whole record.
        keep: When multiple records are found for a lookup, how to return the records.
            - `"first"`: return the first record.
            - `"last"`: return the last record.
            - `False`: return all records.

    Returns:
        A `NamedTuple` of lookup information of the field values with a
        dictionary converter.

    See Also:
        :meth:`~lamindb.models.SQLRecord.search`

    Examples:

        Lookup via auto-complete on `.`::

            import bionty as bt
            bt.Gene.from_source(symbol="ADGB-DT").save()
            lookup = bt.Gene.lookup()
            lookup.adgb_dt

        Look up via auto-complete in dictionary::

            lookup_dict = lookup.dict()
            lookup_dict['ADGB-DT']

        Look up via a specific field::

            lookup_by_ensembl_id = bt.Gene.lookup(field="ensembl_gene_id")
            genes.ensg00000002745

        Return a specific field value instead of the full record::

            lookup_return_symbols = bt.Gene.lookup(field="ensembl_gene_id", return_field="symbol")
    """
    from .sqlrecord import get_name_field

    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    field = get_name_field(registry=queryset.model, field=field)

    return Lookup(
        records=queryset,
        values=[i.get(field) for i in queryset.values()],
        tuple_name=cls.__class__.__name__,
        prefix="ln",
        keep=keep,
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

        Populate the `.parents` ManyToMany relationship (a `QueryManager`)::

            ln.ULabel.from_values(["Label1", "Label2", "Label3"]).save()
            labels = ln.ULabel.filter(name__icontains="label")
            label1 = ln.ULabel.get(name="Label1")
            label1.parents.set(labels)

        Convert all linked parents to a `DataFrame`::

            label1.parents.to_dataframe()
    """

    def to_list(self, field: str | None = None):
        """Populate a list."""
        if field is None:
            return list(self.all())
        else:
            return list(self.values_list(field, flat=True))

    def to_dataframe(self, **kwargs):
        """Convert to DataFrame.

        For `**kwargs`, see :meth:`lamindb.models.QuerySet.to_dataframe`.
        """
        return self.all().to_dataframe(**kwargs)

    @deprecated(new_name="to_dataframe")
    def df(self, **kwargs):
        return self.to_dataframe(**kwargs)

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


# below is just for typing / docs
# Django achieves the same thing with a dynamically generated class
class RelatedManager(QueryManager):
    """Manager for many-to-many and reverse foreign key relationships.

    Provides relationship manipulation methods.

    See Also:
        :class:`lamindb.models.QueryManager`

    Examples:

        Populate the `.parents` ManyToMany relationship (a `RelatedManager`)::

            ln.ULabel.from_values(["Label1", "Label2", "Label3"]).save()
            labels = ln.ULabel.filter(name__icontains="label")
            label1 = ln.ULabel.get(name="Label1")
            label1.parents.set(labels)

        Convert all linked parents to a `DataFrame`::

            label1.parents.to_dataframe()

        Remove a parent label::

            label1.parents.remove(label2)

        Clear all parent labels::

            label1.parents.clear()

    """

    def add(self, *objs, bulk: bool = True) -> None:
        """Add objects to the relationship."""
        ...

    def set(self, objs, *, bulk: bool = True, clear: bool = False) -> None:
        """Set the relationship to the specified objects."""
        ...

    def remove(self, *objs, bulk: bool = True) -> None:
        """Remove objects from the relationship."""
        ...

    def clear(self) -> None:
        """Remove all objects from the relationship."""
        ...
