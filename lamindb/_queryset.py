from datetime import datetime
from typing import Iterable, List, NamedTuple, Optional

import pandas as pd
from django.db import models
from lamindb_setup.dev._docs import doc_args
from lnschema_core import ORM
from lnschema_core.types import ListLike, StrField


class NoResultFound(Exception):
    pass


class MultipleResultsFound(Exception):
    pass


def format_and_convert_to_local_time(series: pd.Series):
    tzinfo = datetime.now().astimezone().tzinfo
    timedelta = tzinfo.utcoffset(datetime.now())  # type: ignore
    return (series + timedelta).dt.strftime("%Y-%m-%d %H:%M:%S")


class QuerySet(models.QuerySet):
    """Extension of Django QuerySet.

    This brings some of the SQLAlchemy/SQLModel/SQL-inspired calls.

    As LaminDB was based on SQLAlchemy/SQLModel in the beginning, and might
    support it again in the future, these calls will be supported longtime.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/4.2/ref/models/querysets/>`__ # noqa

    Examples:

        >>> ln.Tag(name="my tag").save()
        >>> queryset = ln.Tag.select(name="my tag")
        >>> queryset
        <QuerySet [Tag(id=MIeZISeF, name=my tag, updated_at=2023-07-19 19:53:34, created_by_id=DzTjkKse)]> # noqa
    """

    def df(self, include: Optional[List[str]] = None):
        """Convert to ``pd.DataFrame``.

        By default, shows all fields that aren't many-to-many fields, except
        ``created_at``.

        If you'd like to include many-to-many fields, use parameter ``include``.

        Args:
            include: ``Optional[List[str]] = None`` Additional (many-to-many)
                fields to include. Takes expressions like ``"tags__name"``
                ``"cell_types__name"``.

        Examples:

            >>> ln.save(ln.Project.from_values(["Project1", "Project2", "Project3"], field="name")) # noqa
            >>> ln.Project.select().df()
                          name  external_id           updated_at  created_by_id
                  id
            wsCyIq2Z  Project1         None  2023-07-19 19:14:08       DzTjkKse
            MvpDP8Y3  Project2         None  2023-07-19 19:14:08       DzTjkKse
            zKFFabCu  Project3         None  2023-07-19 19:14:08       DzTjkKse
            >>> project = ln.Project.select(name="Project1").one()
            >>> tag = ln.Tag.select(name="benchmark").one()
            >>> project.tags.add(tag)
            >>> ln.Project.select().df(include=["tags__name", "tags__created_by_id"])
                      tags__name  tags__created_by_id      name  external_id           updated_at  created_by_id # noqa
                  id
            wsCyIq2Z  [benchmark]          [DzTjkKse]  Project1         None  2023-07-19 19:14:08       DzTjkKse # noqa
            MvpDP8Y3         None                None  Project2         None  2023-07-19 19:14:08       DzTjkKse # noqa
            zKFFabCu         None                None  Project3         None  2023-07-19 19:14:08       DzTjkKse # noqa
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
        if len(df) > 0 and "updated_at" in df:
            df.updated_at = format_and_convert_to_local_time(df.updated_at)
        if len(df) > 0 and "run_at" in df:
            df.run_at = format_and_convert_to_local_time(df.run_at)
        if "id" in df.columns:
            df = df.set_index("id")
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
                ORM = self.model
                field = getattr(ORM, field_name)
                if not isinstance(field.field, models.ManyToManyField):
                    raise ValueError("Only many-to-many fields are allowed here.")
                related_ORM = (
                    field.field.model
                    if field.field.model != ORM
                    else field.field.related_model
                )
                if field.field.model == related_ORM:
                    left_side_link_model = f"from_{ORM.__name__.lower()}"
                    values_expression = f"to_{ORM.__name__.lower()}__{lookup_str}"
                else:
                    left_side_link_model = f"{ORM.__name__.lower()}"
                    values_expression = f"{related_ORM.__name__.lower()}__{lookup_str}"
                link_df = pd.DataFrame(
                    field.through.objects.values(
                        left_side_link_model, values_expression
                    )
                )
                link_groupby = link_df.groupby(left_side_link_model)[
                    values_expression
                ].apply(list)
                df = pd.concat((link_groupby, df), axis=1)
                df.rename(columns={values_expression: expression}, inplace=True)
        return df

    def list(self, field: Optional[str] = None) -> List[ORM]:
        """Populate a list with the results.

        Examples:

            >>> ln.save(ln.Project.from_values(["Project1", "Project2", "Project3"], field="name")) # noqa
            >>> queryset = ln.Project.select(name__icontains = "project")
            >>> queryset.list()
            [Project(id=NAgTZxoo, name=Project1, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse), # noqa
            Project(id=bnsAgKRC, name=Project2, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse), # noqa
            Project(id=R8xhAJNE, name=Project3, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse)] # noqa
            >>> queryset.list("name")
            ['Project1', 'Project2', 'Project3']
        """
        if field is None:
            return [item for item in self]
        else:
            return [item for item in self.values_list(field, flat=True)]

    def first(self) -> Optional[ORM]:
        """If non-empty, the first result in the query set, otherwise None.

        Examples:
            >>> ln.save(ln.Project.from_values(["Project1", "Project2", "Project3"], field="name")) # noqa
            >>> queryset = ln.Project.select(name__icontains = "project")
            >>> queryset.first()
            Project(id=NAgTZxoo, name=Project1, updated_at=2023-07-19 19:25:48, created_by_id=DzTjkKse) # noqa
        """
        if len(self) == 0:
            return None
        return self[0]

    def one(self) -> ORM:
        """Exactly one result. Throws error if there are more or none.

        Examples:
            >>> ln.Tag(name="benchmark").save()
            >>> ln.Tag.select(name="benchmark").one()
            Tag(id=gznl0GZk, name=benchmark, updated_at=2023-07-19 19:39:01, created_by_id=DzTjkKse) # noqa
        """
        if len(self) == 0:
            raise NoResultFound
        elif len(self) > 1:
            raise MultipleResultsFound
        else:
            return self[0]

    def one_or_none(self) -> Optional[ORM]:
        """At most one result. Returns it if there is one, otherwise returns None.

        Examples:
            >>> ln.Tag(name="benchmark").save()
            >>> ln.Tag.select(name="benchmark").one_or_none()
            Tag(id=gznl0GZk, name=benchmark, updated_at=2023-07-19 19:39:01, created_by_id=DzTjkKse) # noqa
            >>> ln.Tag.select(name="non existing tag").one_or_none()
            None
        """
        if len(self) == 0:
            return None
        elif len(self) == 1:
            return self[0]
        else:
            raise MultipleResultsFound

    @doc_args(ORM.search.__doc__)
    def search(self, string: str, **kwargs):
        """{}"""
        from ._orm import _search

        return _search(cls=self, string=string, **kwargs)

    @doc_args(ORM.lookup.__doc__)
    def lookup(self, field: Optional[StrField] = None) -> NamedTuple:
        """{}"""
        from ._orm import _lookup

        return _lookup(cls=self, field=field)

    @doc_args(ORM.inspect.__doc__)
    def inspect(self, identifiers: ListLike, field: StrField, **kwargs):
        """{}"""
        from ._orm import _inspect

        return _inspect(cls=self, identifiers=identifiers, field=field, **kwargs)

    @doc_args(ORM.map_synonyms.__doc__)
    def map_synonyms(self, synonyms: Iterable, **kwargs):
        """{}"""
        from ._orm import _map_synonyms

        return _map_synonyms(cls=self, synonyms=synonyms, **kwargs)


setattr(models.QuerySet, "df", QuerySet.df)
setattr(models.QuerySet, "list", QuerySet.list)
setattr(models.QuerySet, "first", QuerySet.first)
setattr(models.QuerySet, "one", QuerySet.one)
setattr(models.QuerySet, "one_or_none", QuerySet.one_or_none)
setattr(models.QuerySet, "search", QuerySet.search)
setattr(models.QuerySet, "lookup", QuerySet.lookup)
setattr(models.QuerySet, "inspect", QuerySet.inspect)
setattr(models.QuerySet, "map_synonyms", QuerySet.map_synonyms)
