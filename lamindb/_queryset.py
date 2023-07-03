from typing import List, Optional

from django.db import models


class NoResultFound(Exception):
    pass


class MultipleResultsFound(Exception):
    pass


class QuerySet(models.QuerySet):
    """Extension of Django QuerySet.

    This brings some of the SQLAlchemy/SQLModel/SQL-inspired calls.

    As LaminDB was based on SQLAlchemy/SQLModel in the beginning, and might
    support it again in the future, these calls will be supported longtime.
    """

    def df(self, lists: Optional[List[str]] = None):
        """Return query as DataFrame.

        Include additional fields.
        """
        import pandas as pd

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
            df.updated_at = df.updated_at.dt.strftime("%Y-%m-%d %H:%M:%S")
        if "id" in df.columns:
            df = df.set_index("id")
        if lists is not None:
            if isinstance(lists, str):
                lists = [lists]
            for expression in lists:
                split = expression.split("__")
                field_name = split[0]
                if len(split) > 1:
                    lookup = "__".join(split[1:])
                else:
                    lookup = "id"
                ORM = self.model
                field = getattr(ORM, field_name)
                expression = f"{field.field.model.__name__.lower()}__{lookup}"
                link_df = pd.DataFrame(
                    field.through.objects.values(ORM.__name__.lower(), expression)
                )
                link_groupby = link_df.groupby(ORM.__name__.lower())[expression].apply(
                    list
                )
                df = df.join(link_groupby)
        return df

    def list(self) -> List:
        return list(self)

    def first(self):
        if len(self) == 0:
            return None
        return self[0]

    def one(self):
        if len(self) == 0:
            raise NoResultFound
        elif len(self) > 1:
            raise MultipleResultsFound
        else:
            return self[0]

    def one_or_none(self):
        if len(self) == 0:
            return None
        elif len(self) == 1:
            return self[0]
        else:
            raise MultipleResultsFound
