from __future__ import annotations

from collections import UserList
from typing import TYPE_CHECKING, Iterable, NamedTuple

import pandas as pd
from django.db import models
from django.db.models import F
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import (
    Artifact,
    CanValidate,
    Collection,
    IsVersioned,
    Record,
    Run,
    Transform,
)

from lamindb.core.exceptions import DoesNotExist

if TYPE_CHECKING:
    from lnschema_core.types import ListLike, StrField


class MultipleResultsFound(Exception):
    pass


# def format_and_convert_to_local_time(series: pd.Series):
#     tzinfo = datetime.now().astimezone().tzinfo
#     timedelta = tzinfo.utcoffset(datetime.now())  # type: ignore
#     return (series + timedelta).dt.strftime("%Y-%m-%d %H:%M:%S %Z")


def get_keys_from_df(data: list, registry: Record) -> list[str]:
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
        raise DoesNotExist
    elif len(self) > 1:
        raise MultipleResultsFound(self)
    else:
        return self[0]


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
        return super(models.QuerySet, qs).get(id=idlike)
    elif isinstance(idlike, str):
        qs = qs.filter(uid__startswith=idlike)
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
        return registry.objects.get(**expressions)


class RecordsList(UserList):
    """Is ordered, can't be queried, but has `.df()`."""

    def __init__(self, records: Iterable[Record]):
        super().__init__(record for record in records)

    def df(self) -> pd.DataFrame:
        keys = get_keys_from_df(self.data, self.data[0].__class__)
        values = [record.__dict__ for record in self.data]
        return pd.DataFrame(values, columns=keys)

    def one(self) -> Record:
        """Exactly one result. Throws error if there are more or none."""
        return one_helper(self)


class QuerySet(models.QuerySet, CanValidate):
    """Sets of records returned by queries.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/4.2/ref/models/querysets/>`__ # noqa

    Examples:

        >>> ln.ULabel(name="my label").save()
        >>> queryset = ln.ULabel.filter(name="my label")
        >>> queryset
    """

    @doc_args(Record.df.__doc__)
    def df(
        self, include: str | list[str] | None = None, join: str = "inner"
    ) -> pd.DataFrame:
        """{}"""  # noqa: D415
        # re-order the columns
        exclude_field_names = ["created_at"]
        field_names = [
            field.name
            for field in self.model._meta.fields
            if (
                not isinstance(field, models.ForeignKey)
                and field.name not in exclude_field_names
            )
        ]
        field_names += [
            f"{field.name}_id"
            for field in self.model._meta.fields
            if isinstance(field, models.ForeignKey)
        ]
        for field_name in ["run_id", "created_at", "created_by_id", "updated_at"]:
            if field_name in field_names:
                field_names.remove(field_name)
                field_names.append(field_name)
        if field_names[0] != "uid" and "uid" in field_names:
            field_names.remove("uid")
            field_names.insert(0, "uid")
        # create the dataframe
        df = pd.DataFrame(self.values(), columns=field_names)
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
                Record = self.model
                field = getattr(Record, field_name)
                if isinstance(field.field, models.ManyToManyField):
                    related_ORM = (
                        field.field.model
                        if field.field.model != Record
                        else field.field.related_model
                    )
                    if Record == related_ORM:
                        left_side_link_model = f"from_{Record.__name__.lower()}"
                        values_expression = (
                            f"to_{Record.__name__.lower()}__{lookup_str}"
                        )
                    else:
                        left_side_link_model = f"{Record.__name__.lower()}"
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
                    df = pd.concat((link_groupby, df), axis=1, join=join)
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
                    df = pd.concat((df_anno, df), axis=1, join=join)
        return df

    def delete(self, *args, **kwargs):
        """Delete all records in the query set."""
        # both Transform & Run might reference artifacts
        if self.model in {Artifact, Collection, Transform, Run}:
            for record in self:
                record.delete(*args, **kwargs)
        else:
            self._delete_base_class(*args, **kwargs)

    def list(self, field: str | None = None) -> list[Record]:
        """Populate a list with the results.

        Examples:
            >>> queryset.list()  # list of records
            >>> queryset.list("name")  # list of values
        """
        if field is None:
            return list(self)
        else:
            return list(self.values_list(field, flat=True))

    def first(self) -> Record | None:
        """If non-empty, the first result in the query set, otherwise ``None``.

        Examples:
            >>> queryset.first()
        """
        if len(self) == 0:
            return None
        return self[0]

    def get(self, idlike: int | str | None = None, **expressions) -> Record:
        """Query a single record. Raises error if there are more or none."""
        return get(self, idlike, **expressions)

    def one(self) -> Record:
        """Exactly one result. Raises error if there are more or none."""
        return one_helper(self)

    def one_or_none(self) -> Record | None:
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

    def latest_version(self) -> QuerySet:
        """Filter every version family by latest version."""
        if issubclass(self.model, IsVersioned):
            return self.filter(is_latest=True)
        else:
            raise ValueError("Record isn't subclass of `lamindb.core.IsVersioned`")

    @doc_args(Record.search.__doc__)
    def search(self, string: str, **kwargs):
        """{}"""  # noqa: D415
        from ._record import _search

        return _search(cls=self, string=string, **kwargs)

    @doc_args(Record.lookup.__doc__)
    def lookup(self, field: StrField | None = None, **kwargs) -> NamedTuple:
        """{}"""  # noqa: D415
        from ._record import _lookup

        return _lookup(cls=self, field=field, **kwargs)

    @doc_args(CanValidate.validate.__doc__)
    def validate(self, values: ListLike, field: str | StrField | None = None, **kwargs):
        """{}"""  # noqa: D415
        from ._can_validate import _validate

        return _validate(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanValidate.inspect.__doc__)
    def inspect(self, values: ListLike, field: str | StrField | None = None, **kwargs):
        """{}"""  # noqa: D415
        from ._can_validate import _inspect

        return _inspect(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanValidate.standardize.__doc__)
    def standardize(
        self, values: Iterable, field: str | StrField | None = None, **kwargs
    ):
        """{}"""  # noqa: D415
        from ._can_validate import _standardize

        return _standardize(cls=self, values=values, field=field, **kwargs)


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
