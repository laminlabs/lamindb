from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from django.db.models import QuerySet
from lamindb_setup.dev._docs import doc_args
from lnschema_core.models import Registry, ValidationAware
from lnschema_core.types import ListLike, StrField

from ._from_values import _has_species_field


@classmethod  # type: ignore
@doc_args(ValidationAware.inspect.__doc__)
def inspect(
    cls,
    identifiers: ListLike,
    field: StrField,
    *,
    case_sensitive: bool = False,
    inspect_synonyms: bool = True,
    return_df: bool = False,
    logging: bool = True,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    return _inspect(
        cls=cls,
        identifiers=identifiers,
        field=field,
        case_sensitive=case_sensitive,
        inspect_synonyms=inspect_synonyms,
        return_df=return_df,
        logging=logging,
        **kwargs,
    )


@classmethod  # type: ignore
@doc_args(ValidationAware.validate.__doc__)
def validate(cls, values: ListLike, field: StrField, **kwargs) -> np.ndarray[bool]:
    """{}"""
    if isinstance(values, str):
        values = [values]
    if not isinstance(field, str):
        field = field.field.name

    orm = cls.model if isinstance(cls, QuerySet) else cls
    std_values = set(
        _filter_query_based_on_species(
            orm=orm, species=kwargs.get("species"), values_list_field=field
        )
    )
    return pd.Index(values).isin(std_values)


def _inspect(
    cls,
    identifiers: ListLike,
    field: StrField,
    *,
    case_sensitive: bool = False,
    inspect_synonyms: bool = True,
    return_df: bool = False,
    logging: bool = True,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    from lamin_utils._inspect import inspect

    if not isinstance(field, str):
        field = field.field.name

    orm = cls.model if isinstance(cls, QuerySet) else cls

    return inspect(
        df=_filter_query_based_on_species(orm=orm, species=kwargs.get("species")),
        identifiers=identifiers,
        field=str(field),
        case_sensitive=case_sensitive,
        inspect_synonyms=inspect_synonyms,
        return_df=return_df,
        logging=logging,
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
