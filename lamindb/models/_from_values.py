from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from lamin_utils import colors, logger

if TYPE_CHECKING:
    from lamindb.base.types import FieldAttr, ListLike

    from .query_set import SQLRecordList
    from .sqlrecord import SQLRecord


# The base function for `from_values`
def _from_values(
    iterable: ListLike,
    field: FieldAttr,
    *,
    create: bool = False,
    organism: SQLRecord | str | None = None,
    source: SQLRecord | None = None,
    standardize: bool = True,
    from_source: bool = True,
    mute: bool = False,
    **filter_kwargs,
) -> SQLRecordList:
    """Get or create records from iterables."""
    from .query_set import SQLRecordList

    registry = field.field.model  # type: ignore
    organism_record = get_organism_record_from_field(field, organism, values=iterable)
    # TODO: the create is problematic if field is not a name field
    if create:
        create_kwargs = {}
        if organism_record:
            create_kwargs["organism"] = organism_record
        return SQLRecordList(
            [
                registry(**{field.field.name: value}, **create_kwargs)
                for value in iterable
            ]
        )  # type: ignore

    iterable_idx = index_iterable(iterable)

    # returns existing records & non-existing values
    records, nonexist_values, msg = get_existing_records(
        iterable_idx=iterable_idx,
        field=field,
        organism=organism_record,
        mute=mute,
        **filter_kwargs,
    )

    # new records to be created based on new values
    if len(nonexist_values) > 0:
        if from_source and registry.__base__.__name__ == "BioRecord":
            # if can and needed, get organism record from the existing records
            if (
                organism_record is None
                and len(records) > 0
                and registry.require_organism()
            ):
                organism_record = records[0].organism
            records_public, unmapped_values = create_records_from_source(
                iterable_idx=nonexist_values,
                field=field,
                organism=organism_record,
                source=source,
                standardize=standardize,
                msg=msg,
                mute=mute,
            )
            if len(records_public) > 0:
                msg = ""
            for record in records_public:
                record._from_source = True
            records += records_public
        else:
            unmapped_values = nonexist_values
        # unmapped new_ids will NOT create records
        if len(unmapped_values) > 0:
            # first log the success message
            if len(msg) > 0 and not mute:
                logger.success(msg)
            s = "" if len(unmapped_values) == 1 else "s"
            print_values = colors.yellow(_format_values(unmapped_values))
            n_nonval = colors.yellow(f"{len(unmapped_values)} non-validated")
            if not mute:
                logger.info(
                    f"{colors.red('did not create')} {registry.__name__} record{s} for "
                    f"{n_nonval} {colors.italic(f'{field.field.name}{s}')}: {print_values}"  # type: ignore
                )
    return SQLRecordList(records)


def get_existing_records(
    iterable_idx: pd.Index,
    field: FieldAttr,
    organism: SQLRecord | None = None,
    standardize: bool = True,
    mute: bool = False,
    **filter_kwargs,
) -> tuple[list, pd.Index, str]:
    """Get existing records from the database."""
    from .can_curate import _validate

    # NOTE: existing records matching is agnostic to the source
    registry = field.field.model  # type: ignore
    queryset = registry.filter(**filter_kwargs)

    if standardize:
        # log synonyms mapped terms
        if hasattr(registry, "standardize"):
            syn_mapper = queryset.standardize(
                iterable_idx,
                field=field,
                organism=organism,
                mute=True,
                from_source=False,  # standardize only based on the DB reference
                return_mapper=True,
            )
            iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index
    else:
        syn_mapper = {}

    # now we have to sort the list of queried records
    # preserved = Case(
    #     *[
    #         When(**{field.field.name: value}, then=pos)
    #         for pos, value in enumerate(iterable_idx)
    #     ]
    # )
    # order by causes a factor 10 in runtime
    # records = query_set.order_by(preserved).to_list()

    # log validated terms
    is_validated = _validate(
        cls=queryset, values=iterable_idx, field=field, organism=organism, mute=True
    )
    if len(is_validated) > 0:
        validated = iterable_idx[is_validated]
    else:
        validated = []
    msg = ""
    syn_msg = ""
    if not mute:
        if len(validated) > 0:
            s = "" if len(validated) == 1 else "s"
            print_values = colors.green(_format_values(validated))
            msg = (
                "loaded"
                f" {colors.green(f'{len(validated)} {registry.__name__} record{s}')}"
                f" matching {colors.italic(f'{field.field.name}')}: {print_values}"
            )
        if len(syn_mapper) > 0:
            s = "" if len(syn_mapper) == 1 else "s"
            names = list(syn_mapper.keys())
            print_values = colors.green(_format_values(names))
            syn_msg = (
                "loaded"
                f" {colors.green(f'{len(syn_mapper)} {registry.__name__} record{s}')}"
                f" matching {colors.italic('synonyms')}: {print_values}"
            )

    # no logging if all values are validated
    # logs if there are synonyms
    if len(syn_msg) > 0:
        if len(msg) > 0 and not mute:
            logger.success(msg)
        if not mute:
            logger.success(syn_msg)
        msg = ""

    # get all existing records in the db
    query = {f"{field.field.name}__in": iterable_idx.values}  # type: ignore
    if organism is not None:
        query["organism"] = organism
    records = queryset.filter(**query).to_list()

    if len(validated) == len(iterable_idx):
        return records, pd.Index([]), msg
    else:
        nonval_values = iterable_idx.difference(validated)
        return records, nonval_values, msg


def create_records_from_source(
    iterable_idx: pd.Index,
    field: FieldAttr,
    organism: SQLRecord | None = None,
    source: SQLRecord | None = None,
    standardize: bool = True,
    msg: str = "",
    mute: bool = False,
) -> tuple[list, pd.Index]:
    """Create records from source."""
    registry = field.field.model  # type: ignore
    records: list = []
    # populate additional fields from public_df
    from bionty._organism import OrganismNotSet
    from bionty._source import filter_public_df_columns, get_source_record

    # get the default source
    if organism is None and registry.require_organism(field=field):
        raise OrganismNotSet(
            f"`organism` is required to create new {registry.__name__} records from source!"
        )
    try:
        source_record = get_source_record(registry, organism, source)
    except ValueError:
        # no source found
        return records, iterable_idx

    # create the corresponding PublicOntology object from registry
    try:
        public_ontology = registry.public(source=source_record)
    except Exception:
        # no public source
        return records, iterable_idx

    # filter the columns in public df based on fields
    public_df = filter_public_df_columns(
        registry=registry, public_ontology=public_ontology
    )

    if public_df.empty:
        return records, iterable_idx

    # standardize in the public reference
    # do not inspect synonyms if the field is not name field
    result = public_ontology.inspect(
        iterable_idx,
        field=field.field.name,  # type: ignore
        standardize=False
        if hasattr(registry, "_name_field") and field.field.name != registry._name_field
        else standardize,  # type: ignore
        mute=True,
    )
    syn_mapper = result.synonyms_mapper

    msg_syn: str = ""
    if len(syn_mapper) > 0:
        s = "" if len(syn_mapper) == 1 else "s"
        names = list(syn_mapper.keys())
        print_values = colors.purple(_format_values(names))
        msg_syn = (
            "created"
            f" {colors.purple(f'{len(syn_mapper)} {registry.__name__} record{s} from Bionty')}"
            f" matching {colors.italic('synonyms')}: {print_values}"
        )

        iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # create records for values that are found in the public reference
    # matching either field or synonyms
    mapped_values = iterable_idx.intersection(public_df[field.field.name])  # type: ignore

    multi_msg = ""
    if len(mapped_values) > 0:
        public_kwargs, multi_msg = _bulk_create_dicts_from_df(
            keys=mapped_values,
            column_name=field.field.name,  # type: ignore
            df=public_df,
        )

        create_kwargs = (
            {"organism": organism, "source": source_record}
            if organism is not None
            else {"source": source_record}
        )
        for bk in public_kwargs:
            # skip validation to speed up bulk creation since the values don't validate in the registry DB yet
            records.append(registry(**bk, **create_kwargs, _skip_validation=True))

        # number of records that matches field (not synonyms)
        validated = result.validated
        if len(validated) > 0:
            s = "" if len(validated) == 1 else "s"
            print_values = colors.purple(_format_values(validated))
            # this is the success msg for existing records in the DB from get_existing_records
            if len(msg) > 0 and not mute:
                logger.success(msg)
            if not mute:
                logger.success(
                    "created"
                    f" {colors.purple(f'{len(validated)} {registry.__name__} record{s} from Bionty')}"
                    f" matching {colors.italic(f'{field.field.name}')}: {print_values}"  # type: ignore
                )

    # make sure that synonyms logging appears after the field logging
    if len(msg_syn) > 0 and not mute:
        logger.success(msg_syn)
    # warning about multi matches
    if len(multi_msg) > 0 and not mute:
        logger.warning(multi_msg)

    # return the values that are not found in the public reference
    unmapped_values = iterable_idx.difference(mapped_values)
    return records, unmapped_values


def index_iterable(iterable: ListLike) -> pd.Index:
    """Get unique values from an iterable."""
    idx = pd.Index(iterable).unique()
    # No entries are made for NAs, '', None
    # returns an ordered unique not null list
    return idx[(idx != "") & (~idx.isnull())]


def _format_values(
    names: ListLike, n: int = 20, quotes: bool = True, sep: str = "'"
) -> str:
    """Format values for printing."""
    items = {str(name): None for name in names if name != "None"}

    unique_items = list(items.keys())

    if quotes:
        unique_items = [f"{sep}{item}{sep}" for item in unique_items]

    print_values = ", ".join(unique_items[:n])

    if len(unique_items) > n:
        print_values += ", ..."

    return print_values


def _bulk_create_dicts_from_df(
    keys: set | list, column_name: str, df: pd.DataFrame
) -> tuple[dict, str]:
    """Get fields from a DataFrame for many rows."""
    multi_msg = ""
    if df.index.name != column_name:
        df = df.set_index(column_name).loc[list(keys)]
    if not df.index.is_unique:
        # return all records for multi-matches with a warning
        dup = df.index[df.index.duplicated()].unique().tolist()
        if len(dup) > 0:
            s = "" if len(dup) == 1 else "s"
            print_values = _format_values(dup)
            multi_msg = (
                f"ambiguous validation in Bionty for {len(dup)} record{s}:"
                f" {print_values}"
            )

    return df.reset_index().to_dict(orient="records"), multi_msg


def get_organism_record_from_field(  # type: ignore
    field: FieldAttr,
    organism: str | SQLRecord | None = None,
    values: ListLike = None,
    using_key: str | None = None,
) -> SQLRecord | None:
    """Get organism record based on which field is used in from_values.

    Args:
        field: the field of the registry for from_values
        organism: the organism to get the organism record for
        values: the values passed to from_values
        using_key: the db to get the organism record from

    Returns:
        The organism record if both conditions are met:
            The organism FK is required for the registry
            The field is not unique (e.g. Gene.symbol) or the organism is not None
    """
    registry = field.field.model
    if registry.__base__.__name__ != "BioRecord":
        return None

    from bionty._organism import (
        create_or_get_organism_record,
        infer_organism_from_ensembl_id,
    )

    if values is None:
        values = []

    # if the field is bionty.Gene.ensembl_gene_id, infer organism from ensembl id
    if (
        registry.__get_name_with_module__() == "bionty.Gene"
        and field.field.name == "ensembl_gene_id"
        and len(values) > 0
        and organism is None
    ):
        # Check if values contain bionty.Gene objects with organism field
        from collections.abc import Iterable

        # first check if we have Gene objects
        for v in values:
            # early return to not loop through all values to find a string
            if isinstance(v, str):
                break
            if isinstance(v, registry) and v.organism is not None:
                return v.organism
            # Handle iterables containing Gene objects (but not strings, which are also iterable)
            elif isinstance(v, Iterable) and not isinstance(v, str):
                for item in v:
                    if isinstance(item, registry) and item.organism is not None:
                        return item.organism

        # If no bionty.Gene with organism found, fall back to string-based inference
        # pass the first ensembl id that starts with ENS to infer organism
        first_ensembl = next(
            (v for v in values if isinstance(v, str) and v.startswith("ENS")), ""
        )
        if first_ensembl:
            return infer_organism_from_ensembl_id(first_ensembl, using_key)

    return create_or_get_organism_record(
        organism=organism, registry=registry, field=field
    )
