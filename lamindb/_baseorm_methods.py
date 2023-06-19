import builtins
from typing import NamedTuple, Optional, Union

from django.db.models import CharField, TextField
from lamin_logger import logger
from lamin_logger._lookup import Lookup
from lnschema_core import BaseORM
from pandas import DataFrame

_is_ipython = getattr(builtins, "__IPYTHON__", False)


class ValidationError(Exception):
    pass


def validate_required_fields(orm: BaseORM, kwargs):
    required_fields = {
        k.name for k in orm._meta.fields if not k.null and k.default is None
    }
    required_fields_not_passed = {k: None for k in required_fields if k not in kwargs}
    kwargs.update(required_fields_not_passed)
    missing_fields = [
        k for k, v in kwargs.items() if v is None and k in required_fields
    ]
    if missing_fields:
        raise TypeError(f"{missing_fields} are required.")


def suggest_objects_with_same_name(orm: BaseORM, kwargs) -> Optional[str]:
    if "name" not in kwargs:
        return None
    elif kwargs["name"] is None:
        return None
    else:
        try:
            results = orm.search(kwargs["name"])
            # let's subset results to those with at least 0.5 levensteihn distance
            results = results.loc[results.__ratio__ >= 0.5]
        except KeyError:  # will be fixed soon
            return None
        # test for exact match
        if len(results) > 0:
            if results.index[0] == kwargs["name"]:
                logger.warning("Object with exact same name exists, returning it")
                return "object-with-same-name-exists"
            else:
                msg = "Entries with similar names exist:"
                if _is_ipython:
                    from IPython.display import display

                    logger.warning(f"{msg}")
                    display(results)
                else:
                    logger.warning(f"{msg}\n{results.name}")
    return None


def __init__(orm: BaseORM, *args, **kwargs):
    if not args:  # if args, object is loaded from DB
        validate_required_fields(orm, kwargs)
        result = suggest_objects_with_same_name(orm, kwargs)
        if result == "object-with-same-name-exists":
            existing_object = orm.select(name=kwargs["name"])[0]
            new_args = [
                getattr(existing_object, field.attname)
                for field in orm._meta.concrete_fields
            ]
            super(BaseORM, orm).__init__(*new_args)
            orm._state.adding = False  # mimic from_db
            return None
    super(BaseORM, orm).__init__(*args, **kwargs)


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
        field: `Optional[Union[str, CharField, TextField]] = None` The field to
            look up the values for. Defaults to 'name'.

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


BaseORM.__init__ = __init__
BaseORM.search = search
BaseORM.lookup = lookup
