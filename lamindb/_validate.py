from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from django.db.models import QuerySet
from lamin_utils import colors, logger
from lamin_utils._inspect import InspectResult
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
    mute: bool = False,
    **kwargs,
) -> InspectResult:
    """{}"""
    return _inspect(
        cls=cls,
        values=values,
        field=field,
        mute=mute,
        **kwargs,
    )


@classmethod  # type: ignore
@doc_args(ValidationAware.validate.__doc__)
def validate(
    cls, values: ListLike, field: StrField, *, mute: bool = False, **kwargs
) -> np.ndarray:
    """{}"""
    return _validate(cls=cls, values=values, field=field, mute=mute, **kwargs)


def _inspect(
    cls,
    values: ListLike,
    field: StrField,
    *,
    mute: bool = False,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    from lamin_utils._inspect import inspect

    if isinstance(values, str):
        values = [values]
    if not isinstance(field, str):
        field = field.field.name

    orm = cls.model if isinstance(cls, QuerySet) else cls

    # inspect in the DB
    result_db = inspect(
        df=_filter_query_based_on_species(orm=orm, species=kwargs.get("species")),
        identifiers=values,
        field=str(field),
        mute=mute,
        **kwargs,
    )
    nonval = set(result_db.non_validated).difference(result_db.synonyms_mapper.keys())

    if len(nonval) > 0 and orm.__get_schema_name__() == "bionty":
        bionty_result = orm.bionty().inspect(
            values=nonval, field=str(field), mute=True, **kwargs
        )
        bionty_validated = bionty_result.validated
        bionty_mapper = bionty_result.synonyms_mapper
        hint = False
        if len(bionty_validated) > 0 and not mute:
            print_values = ", ".join(bionty_validated[:20])
            if len(bionty_validated) > 20:
                print_values += ", ..."
            s = "" if len(bionty_validated) == 1 else "s"
            logger.info(
                f"-- detected {len(bionty_validated)} term{s} in Bionty for"
                f" {str(field)}: {print_values}"
            )
            hint = True

        if len(bionty_mapper) > 0 and not mute:
            print_values = ", ".join(list(bionty_mapper.keys())[:20])
            if len(bionty_mapper) > 20:
                print_values += ", ..."
            s = "" if len(bionty_mapper) == 1 else "s"
            logger.info(
                f"-- detected {len(bionty_mapper)} term{s} in Bionty as synonym{s}:"
                f" {print_values}"
            )
            hint = True

        if hint:
            logger.hint(
                "   add records from Bionty to your registry via"
                f" {colors.italic('.from_values()')}"
            )

    return result_db


def _validate(
    cls, values: ListLike, field: StrField, *, mute: bool = False, **kwargs
) -> np.ndarray:
    """{}"""
    from lamin_utils._inspect import validate

    if isinstance(values, str):
        values = [values]
    if not isinstance(field, str):
        field = field.field.name

    orm = cls.model if isinstance(cls, QuerySet) else cls
    field_values = pd.Series(
        _filter_query_based_on_species(
            orm=orm,
            species=kwargs.get("species"),
            values_list_field=field,
        ),
        dtype="object",
    )
    return validate(
        identifiers=values,
        field_values=field_values,
        case_sensitive=True,
        mute=mute,
        field=field,
        **kwargs,
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
