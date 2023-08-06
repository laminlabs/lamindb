from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from django.db.models import QuerySet
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Registry, ValidationAware
from lnschema_core.types import ListLike, StrField

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import _has_species_field


@classmethod  # type: ignore
@doc_args(ValidationAware.inspect.__doc__)
def inspect(
    cls,
    values: ListLike,
    field: StrField,
    *,
    return_df: bool = False,
    mute: bool = False,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    return _inspect(
        cls=cls,
        values=values,
        field=field,
        return_df=return_df,
        mute=mute,
        **kwargs,
    )


@classmethod  # type: ignore
@doc_args(ValidationAware.validate.__doc__)
def validate(cls, values: ListLike, field: StrField, **kwargs) -> np.ndarray[bool]:
    """{}"""
    from lamin_utils._inspect import validate

    if isinstance(values, str):
        values = [values]
    if not isinstance(field, str):
        field = field.field.name

    orm = cls.model if isinstance(cls, QuerySet) else cls
    field_values = pd.Series(
        _filter_query_based_on_species(
            orm=orm, species=kwargs.get("species"), values_list_field=field
        )
    )
    return validate(
        values=values,
        field_values=field_values,
        case_sensitive=True,
        return_df=False,
    )


def _inspect(
    cls,
    values: ListLike,
    field: StrField,
    *,
    return_df: bool = False,
    mute: bool = False,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    from lamin_utils._inspect import inspect

    if not isinstance(field, str):
        field = field.field.name

    orm = cls.model if isinstance(cls, QuerySet) else cls

    return inspect(
        df=_filter_query_based_on_species(orm=orm, species=kwargs.get("species")),
        values=values,
        field=str(field),
        inspect_synonyms=True,
        return_df=return_df,
        logging=not mute,
    )


def _filter_query_based_on_species(
    orm: Union[Registry, QuerySet],
    species: Optional[Union[str, Registry]] = None,
    values_list_field: Optional[str] = None,
):
    import pandas as pd

    if values_list_field is None:
        records = orm.all() if isinstance(orm, QuerySet) else orm.objects.all()
    else:
        records = orm if isinstance(orm, QuerySet) else orm.objects
    if _has_species_field(orm):
        # here, we can safely import lnschema_bionty
        from lnschema_bionty._bionty import create_or_get_species_record

        species_record = create_or_get_species_record(
            species=species, orm=orm.model if isinstance(orm, QuerySet) else orm
        )
        if species_record is not None:
            records = records.filter(species__name=species_record.name)

    if values_list_field is None:
        return pd.DataFrame.from_records(records.values())
    else:
        return records.values_list(values_list_field, flat=True)


METHOD_NAMES = [
    "validate",
    "inspect",
]

if _TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(ValidationAware, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, ValidationAware, globals())
