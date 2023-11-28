from typing import Iterable, List, NamedTuple, Optional, Union

import pandas as pd
from django.db import models
from lamindb_setup.dev._docs import doc_args
from lnschema_core.models import CanValidate, IsTree, Registry
from lnschema_core.types import ListLike, StrField


class NoResultFound(Exception):
    pass


class MultipleResultsFound(Exception):
    pass


# def format_and_convert_to_local_time(series: pd.Series):
#     tzinfo = datetime.now().astimezone().tzinfo
#     timedelta = tzinfo.utcoffset(datetime.now())  # type: ignore
#     return (series + timedelta).dt.strftime("%Y-%m-%d %H:%M:%S %Z")


class QuerySet(models.QuerySet, CanValidate, IsTree):
    """Lazily loaded queried records returned by queries.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/4.2/ref/models/querysets/>`__ # noqa

    Examples:

        >>> ln.ULabel(name="my label").save()
        >>> queryset = ln.ULabel.filter(name="my label")
        >>> queryset
        <QuerySet [ULabel(id=MIeZISeF, name=my label, updated_at=2023-07-19 19:53:34, created_by_id=DzTjkKse)]> # noqa
    """

    def df(self, include: Optional[List[str]] = None) -> pd.DataFrame:
        """Convert to ``pd.DataFrame``.

        By default, shows all fields that aren't many-to-many fields, except
        ``created_at``.

        If you'd like to include many-to-many fields, use parameter ``include``.

        Args:
            include: ``Optional[List[str]] = None`` Additional (many-to-many)
                fields to include. Takes expressions like ``"labels__name"``
                ``"cell_types__name"``.

        Examples:

            >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name")) # noqa
            >>> ln.ULabel.filter().df()
            >>> label = ln.ULabel.filter(name="ULabel1").one()
            >>> label = ln.ULabel.filter(name="benchmark").one()
            >>> label.parents.add(label)
            >>> ln.ULabel.filter().df(include=["labels__name", "labels__created_by_id"])
        """
        data = self.values()
        if len(data) > 0:
            keys = list(data[0].keys())
            if "created_at" in keys:
                keys.remove("created_at")
        else:
            keys = [
                field.name
                for field in self.model._meta.fields
                if (
                    not isinstance(field, models.ForeignKey)
                    and field.name != "created_at"
                )
            ]
            keys += [
                f"{field.name}_id"
                for field in self.model._meta.fields
                if isinstance(field, models.ForeignKey)
            ]
        df = pd.DataFrame(self.values(), columns=keys)
        # if len(df) > 0 and "updated_at" in df:
        #     df.updated_at = format_and_convert_to_local_time(df.updated_at)
        # if len(df) > 0 and "run_at" in df:
        #     df.run_at = format_and_convert_to_local_time(df.run_at)
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
                if not isinstance(field.field, models.ManyToManyField):
                    raise ValueError("Only many-to-many fields are allowed here.")
                related_ORM = (
                    field.field.model
                    if field.field.model != Registry
                    else field.field.related_model
                )
                if Registry == related_ORM:
                    left_side_link_model = f"from_{Registry.__name__.lower()}"
                    values_expression = f"to_{Registry.__name__.lower()}__{lookup_str}"
                else:
                    left_side_link_model = f"{Registry.__name__.lower()}"
                    values_expression = f"{related_ORM.__name__.lower()}__{lookup_str}"
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
        return df

    def list(self, field: Optional[str] = None) -> List[Registry]:
        """Populate a list with the results.

        Examples:

            >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name")) # noqa
            >>> queryset = ln.ULabel.filter(name__icontains = "project")
            >>> queryset.list()
            [ULabel(id=NAgTZxoo, name=ULabel1, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse), # noqa
            ULabel(id=bnsAgKRC, name=ULabel2, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse), # noqa
            ULabel(id=R8xhAJNE, name=ULabel3, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse)] # noqa
            >>> queryset.list("name")
            ['ULabel1', 'ULabel2', 'ULabel3']
        """
        if field is None:
            return [item for item in self]
        else:
            return [item for item in self.values_list(field, flat=True)]

    def first(self) -> Optional[Registry]:
        """If non-empty, the first result in the query set, otherwise None.

        Examples:
            >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name")) # noqa
            >>> queryset = ln.ULabel.filter(name__icontains = "project")
            >>> queryset.first()
            ULabel(id=NAgTZxoo, name=ULabel1, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse) # noqa
        """
        if len(self) == 0:
            return None
        return self[0]

    def one(self) -> Registry:
        """Exactly one result. Throws error if there are more or none.

        Examples:
            >>> ln.ULabel(name="benchmark").save()
            >>> ln.ULabel.filter(name="benchmark").one()
            ULabel(id=gznl0GZk, name=benchmark, updated_at=2023-07-19 19:39:01, created_by_id=DzTjkKse) # noqa
        """
        if len(self) == 0:
            raise NoResultFound
        elif len(self) > 1:
            raise MultipleResultsFound
        else:
            return self[0]

    def one_or_none(self) -> Optional[Registry]:
        """At most one result. Returns it if there is one, otherwise returns None.

        Examples:
            >>> ln.ULabel(name="benchmark").save()
            >>> ln.ULabel.filter(name="benchmark").one_or_none()
            ULabel(id=gznl0GZk, name=benchmark, updated_at=2023-07-19 19:39:01, created_by_id=DzTjkKse) # noqa
            >>> ln.ULabel.filter(name="non existing label").one_or_none()
            None
        """
        if len(self) == 0:
            return None
        elif len(self) == 1:
            return self[0]
        else:
            raise MultipleResultsFound

    @doc_args(Registry.search.__doc__)
    def search(self, string: str, **kwargs):
        """{}"""
        from ._registry import _search

        return _search(cls=self, string=string, **kwargs)

    @doc_args(Registry.lookup.__doc__)
    def lookup(self, field: Optional[StrField] = None, **kwargs) -> NamedTuple:
        """{}"""
        from ._registry import _lookup

        return _lookup(cls=self, field=field, **kwargs)

    @doc_args(CanValidate.validate.__doc__)
    def validate(
        self, values: ListLike, field: Optional[Union[str, StrField]] = None, **kwargs
    ):
        """{}"""
        from ._validate import _validate

        return _validate(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanValidate.inspect.__doc__)
    def inspect(
        self, values: ListLike, field: Optional[Union[str, StrField]] = None, **kwargs
    ):
        """{}"""
        from ._validate import _inspect

        return _inspect(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanValidate.standardize.__doc__)
    def standardize(
        self, values: Iterable, field: Optional[Union[str, StrField]] = None, **kwargs
    ):
        """{}"""
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
        """{}"""
        from .dev._view_tree import view_tree as _view_tree

        _view_tree(
            cls=self,
            level=level,
            limit_to_directories=limit_to_directories,
            length_limit=length_limit,
            max_files_per_dir_per_type=max_files_per_dir_per_type,
        )


setattr(models.QuerySet, "df", QuerySet.df)
setattr(models.QuerySet, "list", QuerySet.list)
setattr(models.QuerySet, "first", QuerySet.first)
setattr(models.QuerySet, "one", QuerySet.one)
setattr(models.QuerySet, "one_or_none", QuerySet.one_or_none)
setattr(models.QuerySet, "search", QuerySet.search)
setattr(models.QuerySet, "lookup", QuerySet.lookup)
setattr(models.QuerySet, "validate", QuerySet.validate)
setattr(models.QuerySet, "inspect", QuerySet.inspect)
setattr(models.QuerySet, "standardize", QuerySet.standardize)
