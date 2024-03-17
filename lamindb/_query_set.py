from collections import UserList
from typing import Dict, Iterable, List, NamedTuple, Optional, Union

import pandas as pd
from django.db import models
from django.db.models import F
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import (
    Artifact,
    CanValidate,
    Collection,
    IsTree,
    IsVersioned,
    Registry,
    Transform,
)
from lnschema_core.types import ListLike, StrField


class NoResultFound(Exception):
    pass


class MultipleResultsFound(Exception):
    pass


# def format_and_convert_to_local_time(series: pd.Series):
#     tzinfo = datetime.now().astimezone().tzinfo
#     timedelta = tzinfo.utcoffset(datetime.now())  # type: ignore
#     return (series + timedelta).dt.strftime("%Y-%m-%d %H:%M:%S %Z")


def get_keys_from_df(data: List, registry: Registry) -> List[str]:
    if len(data) > 0:
        if isinstance(data[0], dict):
            keys = list(data[0].keys())
        else:
            keys = list(data[0].__dict__.keys())
            if "_state" in keys:
                keys.remove("_state")
    else:
        keys = [
            field.name
            for field in registry._meta.fields
            if not isinstance(field, models.ForeignKey)
        ]
        keys += [
            f"{field.name}_id"
            for field in registry._meta.fields
            if isinstance(field, models.ForeignKey)
        ]
    return keys


def one_helper(self):
    if len(self) == 0:
        raise NoResultFound
    elif len(self) > 1:
        raise MultipleResultsFound(self)
    else:
        return self[0]


class RecordsList(UserList):
    """Is ordered, can't be queried, but has `.df()`."""

    def __init__(self, records: Iterable[Registry]):
        super().__init__(record for record in records)

    def df(self) -> pd.DataFrame:
        keys = get_keys_from_df(self.data, self.data[0].__class__)
        values = [record.__dict__ for record in self.data]
        return pd.DataFrame(values, columns=keys)

    def one(self) -> Registry:
        """Exactly one result. Throws error if there are more or none."""
        return one_helper(self)


class QuerySet(models.QuerySet, CanValidate, IsTree):
    """Sets of records returned by queries.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/4.2/ref/models/querysets/>`__ # noqa

    Examples:

        >>> ln.ULabel(name="my label").save()
        >>> queryset = ln.ULabel.filter(name="my label")
        >>> queryset
    """

    @doc_args(Registry.df.__doc__)
    def df(self, include: Optional[Union[str, List[str]]] = None) -> pd.DataFrame:
        """{}."""
        data = self.values()
        keys = get_keys_from_df(data, self.model)
        df = pd.DataFrame(self.values(), columns=keys)
        # if len(df) > 0 and "updated_at" in df:
        #     df.updated_at = format_and_convert_to_local_time(df.updated_at)
        # if len(df) > 0 and "started_at" in df:
        #     df.started_at = format_and_convert_to_local_time(df.started_at)
        pk_name = self.model._meta.pk.name
        pk_column_name = pk_name if pk_name in df.columns else f"{pk_name}_id"
        if pk_column_name in df.columns:
            df = df.set_index(pk_column_name)
        if len(df) == 0:
            return df
        if include is not None:
            if isinstance(include, str):
                include = [include]
            # fix ordering
            include = include[::-1]
            for expression in include:
                split = expression.split("__")
                field_name = split[0]
                if len(split) > 1:
                    lookup_str = "__".join(split[1:])
                else:
                    lookup_str = "id"
                Registry = self.model
                field = getattr(Registry, field_name)
                if isinstance(field.field, models.ManyToManyField):
                    related_ORM = (
                        field.field.model
                        if field.field.model != Registry
                        else field.field.related_model
                    )
                    if Registry == related_ORM:
                        left_side_link_model = f"from_{Registry.__name__.lower()}"
                        values_expression = (
                            f"to_{Registry.__name__.lower()}__{lookup_str}"
                        )
                    else:
                        left_side_link_model = f"{Registry.__name__.lower()}"
                        values_expression = (
                            f"{related_ORM.__name__.lower()}__{lookup_str}"
                        )
                    link_df = pd.DataFrame(
                        field.through.objects.values(
                            left_side_link_model, values_expression
                        )
                    )
                    if link_df.shape[0] == 0:
                        return df
                    link_groupby = link_df.groupby(left_side_link_model)[
                        values_expression
                    ].apply(list)
                    df = pd.concat((link_groupby, df), axis=1, join="inner")
                    df.rename(columns={values_expression: expression}, inplace=True)
                else:
                    # the F() based implementation could also work for many-to-many,
                    # would need to test what is faster
                    df_anno = pd.DataFrame(
                        self.annotate(expression=F(expression)).values(
                            pk_column_name, "expression"
                        )
                    )
                    df_anno = df_anno.set_index(pk_column_name)
                    df_anno.rename(columns={"expression": expression}, inplace=True)
                    df = pd.concat((df_anno, df), axis=1, join="inner")
        return df

    def delete(self, *args, **kwargs):
        """Delete all records in the query set."""
        if self.model in {Artifact, Collection, Transform}:
            for record in self:
                record.delete(*args, **kwargs)
        else:
            self._delete_base_class(*args, **kwargs)

    def list(self, field: Optional[str] = None) -> List[Registry]:
        """Populate a list with the results.

        Examples:
            >>> queryset.list()  # list of records
            >>> queryset.list("name")  # list of values
        """
        if field is None:
            return list(self)
        else:
            return list(self.values_list(field, flat=True))

    def first(self) -> Optional[Registry]:
        """If non-empty, the first result in the query set, otherwise ``None``.

        Examples:
            >>> queryset.first()
        """
        if len(self) == 0:
            return None
        return self[0]

    def one(self) -> Registry:
        """Exactly one result. Raises error if there are more or none.

        Examples:
            >>> ln.ULabel.filter(name="benchmark").one()
        """
        return one_helper(self)

    def one_or_none(self) -> Optional[Registry]:
        """At most one result. Returns it if there is one, otherwise returns ``None``.

        Examples:
            >>> ln.ULabel.filter(name="benchmark").one_or_none()
            >>> ln.ULabel.filter(name="non existing label").one_or_none()
        """
        if len(self) == 0:
            return None
        elif len(self) == 1:
            return self[0]
        else:
            raise MultipleResultsFound(self.all())

    def latest_version(self) -> RecordsList:
        """Filter every version family by latest version."""
        if issubclass(self.model, IsVersioned):
            return filter_query_set_by_latest_version(self)
        else:
            raise ValueError("Registry isn't subclass of `lamindb.core.IsVersioned`")

    @doc_args(Registry.search.__doc__)
    def search(self, string: str, **kwargs):
        """{}."""
        from ._registry import _search

        return _search(cls=self, string=string, **kwargs)

    @doc_args(Registry.lookup.__doc__)
    def lookup(self, field: Optional[StrField] = None, **kwargs) -> NamedTuple:
        """{}."""
        from ._registry import _lookup

        return _lookup(cls=self, field=field, **kwargs)

    @doc_args(CanValidate.validate.__doc__)
    def validate(
        self, values: ListLike, field: Optional[Union[str, StrField]] = None, **kwargs
    ):
        """{}."""
        from ._validate import _validate

        return _validate(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanValidate.inspect.__doc__)
    def inspect(
        self, values: ListLike, field: Optional[Union[str, StrField]] = None, **kwargs
    ):
        """{}."""
        from ._validate import _inspect

        return _inspect(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanValidate.standardize.__doc__)
    def standardize(
        self, values: Iterable, field: Optional[Union[str, StrField]] = None, **kwargs
    ):
        """{}."""
        from ._validate import _standardize

        return _standardize(cls=self, values=values, field=field, **kwargs)

    @doc_args(IsTree.view_tree.__doc__)
    def view_tree(
        self,
        level: int = -1,
        limit_to_directories: bool = False,
        length_limit: int = 1000,
        max_files_per_dir_per_type: int = 7,
    ) -> None:
        """{}."""
        from .core._view_tree import view_tree as _view_tree

        _view_tree(
            cls=self,
            level=level,
            limit_to_directories=limit_to_directories,
            length_limit=length_limit,
            max_files_per_dir_per_type=max_files_per_dir_per_type,
        )


def filter_query_set_by_latest_version(ordered_query_set: QuerySet) -> RecordsList:
    if len(ordered_query_set) == 0:
        return ordered_query_set
    first_record = ordered_query_set[0]
    records_in_view = {}
    records_in_view[first_record.stem_uid] = first_record
    for record in ordered_query_set:
        # this overwrites user-provided ordering (relevant records ordered by a
        # certain field will not show if they are not the latest version)
        if record.stem_uid not in records_in_view:
            records_in_view[record.stem_uid] = record
        else:
            if record.created_at > records_in_view[record.stem_uid].created_at:
                # deleting the entry is needed to preserve the integrity of
                # user-provided ordering
                del records_in_view[record.stem_uid]
                records_in_view[record.stem_uid] = record
    list_records_in_view = RecordsList(records_in_view.values())
    return list_records_in_view


models.QuerySet.df = QuerySet.df
models.QuerySet.list = QuerySet.list
models.QuerySet.first = QuerySet.first
models.QuerySet.one = QuerySet.one
models.QuerySet.one_or_none = QuerySet.one_or_none
models.QuerySet.latest_version = QuerySet.latest_version
models.QuerySet.search = QuerySet.search
models.QuerySet.lookup = QuerySet.lookup
models.QuerySet.validate = QuerySet.validate
models.QuerySet.inspect = QuerySet.inspect
models.QuerySet.standardize = QuerySet.standardize
models.QuerySet._delete_base_class = models.QuerySet.delete
models.QuerySet.delete = QuerySet.delete
