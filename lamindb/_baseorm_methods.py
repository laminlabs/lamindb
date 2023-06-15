from typing import Any, NamedTuple, Union

from django.db import models
from lnschema_core import BaseORM


@classmethod  # type: ignore
def search(
    cls,
    string: str,
    field: Union[models.CharField, str, models.TextField, None] = None,
    synonyms_field: Union[models.TextField, models.CharField, str, None] = "synonyms",
    case_sensitive: bool = True,
    return_ranked_results: bool = False,
    synonyms_sep: str = "|",
) -> Any:
    """Search a given string against a field.

    Args:
        string: The input string to match against the field ontology values.
        field: The field against which the input string is matching.
        synonyms_field: Also map against in the synonyms
            If None, no mapping against synonyms.
        case_sensitive: Whether the match is case sensitive.
        return_ranked_results: If True, return all entries ranked by matching ratios.

    Returns:
        Best match record of the input string.
    """
    import pandas as pd
    from lamin_logger._search import search

    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.objects
    df = pd.DataFrame.from_records(records.values())

    res = search(
        df=df,
        string=string,
        field=field,
        synonyms_field=str(synonyms_field),
        case_sensitive=case_sensitive,
        return_ranked_results=return_ranked_results,
        synonyms_sep=synonyms_sep,
        tuple_name=cls.__name__,
    )

    if return_ranked_results is True or res is None:
        return res
    else:
        if isinstance(res, list):
            return [records.get(id=r.id) for r in res]
        else:
            return records.get(id=res.id)


@classmethod  # type: ignore
def lookup(
    cls, field: Union[models.CharField, models.TextField, str, None] = None
) -> NamedTuple:
    from lamin_logger._lookup import Lookup

    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.objects

    return Lookup(
        records=records.all(),
        values=[i.get(field) for i in records.values()],
        tuple_name=cls.__name__,
        prefix="ln",
    ).lookup()


def get_default_str_field(orm: BaseORM) -> str:
    model_field_names = [i.name for i in orm._meta.fields]

    # set default field
    if "name" in model_field_names:
        # by default use the name field
        field = orm._meta.get_field("name")
    else:
        # first char or text field that doesn't contain "id"
        for i in orm._meta.fields:
            if "id" in i.name:
                continue
            if i.get_internal_type() in {"CharField", "TextField"}:
                field = i
                break

    # no default field can be found
    if field is None:
        raise ValueError("Please specify a field to search against!")

    return field.name


BaseORM.search = search
BaseORM.lookup = lookup
