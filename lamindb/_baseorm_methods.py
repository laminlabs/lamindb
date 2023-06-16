from typing import NamedTuple, Optional, Union

from django.db.models import CharField, TextField
from lamin_logger._lookup import Lookup
from lnschema_core import BaseORM
from pandas import DataFrame


@classmethod  # type: ignore
def search(
    cls,
    string: str,
    *,
    field: Optional[Union[str, CharField, TextField]] = None,
    top_hit: bool = False,
    case_sensitive: bool = True,
    synonyms_field: Optional[Union[str, TextField, CharField]] = "synonyms",
    synonyms_sep: str = "|",
) -> Union[DataFrame, BaseORM]:
    """Search the table.

    Args:
        string: `str` The input string to match against the field ontology values.
        field: `Optional[Union[str, CharField, TextField]] = None` The field
            against which the input string is matching.
        top_hit: `bool = False` If `True`, return only the top hit or hits (in
            case of equal scores).
        case_sensitive: `bool = False` Whether the match is case sensitive.
        synonyms_field: `bool = True` Also search synonyms. If `None`, is ignored.

    Returns:
        Best match record of the input string.
    """
    import pandas as pd
    from lamin_logger._search import search

    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.objects.all()
    df = pd.DataFrame.from_records(records.values())

    result = search(
        df=df,
        string=string,
        field=field,
        synonyms_field=str(synonyms_field),
        case_sensitive=case_sensitive,
        return_ranked_results=not top_hit,
        synonyms_sep=synonyms_sep,
        tuple_name=cls.__name__,
    )

    if not top_hit or result is None:
        return result
    else:
        if isinstance(result, list):
            return [records.get(id=r.id) for r in result]
        else:
            return records.get(id=result.id)


@classmethod  # type: ignore
def lookup(cls, field: Optional[Union[str, CharField, TextField]] = None) -> NamedTuple:
    """Return an auto-complete object for a field.

    Args:
        field: The field to look up the values for.
            Defaults to 'name'.

    Returns:
        A `NamedTuple` of lookup information of the field values with a
        dictionary converter.

    Examples:
        >>> import lnschema_bionty as lb
        >>> lookup = lb.Gene.lookup()
        >>> lookup.adgb_dt
        >>> lookup_dict = lookup.dict()
        >>> lookup['ADGB-DT']
    """
    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.objects.all()

    return Lookup(
        records=records,
        values=[i.get(field) for i in records.values()],
        tuple_name=cls.__name__,
        prefix="ln",
    ).lookup()


lookup.__doc__ = Lookup.__doc__


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
