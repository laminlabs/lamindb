from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from lamin_utils import colors, logger

if TYPE_CHECKING:
    from collections.abc import Iterable

    from lamindb.base.types import FieldAttr, ListLike, StrField

    from .query_set import RecordList
    from .record import Record


# The base function for `from_values`
def get_or_create_records(
    iterable: ListLike,
    field: StrField,
    *,
    create: bool = False,
    organism: Record | str | None = None,
    source: Record | None = None,
    mute: bool = False,
    _from_source: bool = False,
) -> RecordList:
    """Get or create records from iterables."""
    from .query_set import RecordList

    registry = field.field.model  # type: ignore
    if isinstance(field, str):
        field = registry._meta.get_field(field)
    organism_record = _get_organism_record(field, organism, values=iterable)
    # TODO: the create is problematic if field is not a name field
    if create:
        create_kwargs = {}
        if organism_record:
            create_kwargs["organism"] = organism_record
        return RecordList(
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
    )

    # new records to be created based on new values
    if len(nonexist_values) > 0:
        if _from_source and hasattr(registry, "source_id"):
            # if can and needed, get organism record from the existing records
            if (
                organism_record is None
                and len(records) > 0
                and _require_organism(registry)
            ):
                organism_record = records[0].organism
            records_bionty, unmapped_values = create_records_from_source(
                iterable_idx=nonexist_values,
                field=field,
                organism=organism_record,
                source=source,
                msg=msg,
                mute=mute,
            )
            if len(records_bionty) > 0:
                msg = ""
            for record in records_bionty:
                record._from_source = True
            records += records_bionty
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
                logger.warning(
                    f"{colors.red('did not create')} {registry.__name__} record{s} for "
                    f"{n_nonval} {colors.italic(f'{field.field.name}{s}')}: {print_values}"  # type: ignore
                )
    return RecordList(records)


def get_existing_records(
    iterable_idx: pd.Index,
    field: FieldAttr,
    organism: Record | None = None,
    mute: bool = False,
) -> tuple[list, pd.Index, str]:
    # NOTE: existing records matching is agnostic to the source
    model = field.field.model  # type: ignore

    # log synonyms mapped terms
    syn_mapper = model.standardize(
        iterable_idx,
        field=field,
        organism=organism,
        mute=True,
        public_aware=False,  # standardize only based on the DB reference
        return_mapper=True,
    )
    iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # now we have to sort the list of queried records
    # preserved = Case(
    #     *[
    #         When(**{field.field.name: value}, then=pos)
    #         for pos, value in enumerate(iterable_idx)
    #     ]
    # )
    # order by causes a factor 10 in runtime
    # records = query_set.order_by(preserved).list()

    # log validated terms
    is_validated = model.validate(
        iterable_idx, field=field, organism=organism, mute=True
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
                f" {colors.green(f'{len(validated)} {model.__name__} record{s}')}"
                f" matching {colors.italic(f'{field.field.name}')}: {print_values}"
            )
        if len(syn_mapper) > 0:
            s = "" if len(syn_mapper) == 1 else "s"
            names = list(syn_mapper.keys())
            print_values = colors.green(_format_values(names))
            syn_msg = (
                "loaded"
                f" {colors.green(f'{len(syn_mapper)} {model.__name__} record{s}')}"
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
    records = model.filter(**query).list()

    if len(validated) == len(iterable_idx):
        return records, pd.Index([]), msg
    else:
        nonval_values = iterable_idx.difference(validated)
        return records, nonval_values, msg


def create_records_from_source(
    iterable_idx: pd.Index,
    field: FieldAttr,
    organism: Record | None = None,
    source: Record | None = None,
    msg: str = "",
    mute: bool = False,
) -> tuple[list, pd.Index]:
    model = field.field.model  # type: ignore
    records: list = []
    # populate additional fields from bionty
    from bionty._bionty import get_source_record
    from bionty.core._bionty import filter_bionty_df_columns

    # create the corresponding bionty object from model
    try:
        # TODO: more generic
        public_ontology = model.public(organism=organism, source=source)
    except Exception:
        # for custom records that are not created from public sources
        return records, iterable_idx
    # get the default source
    if source is None:
        source = get_source_record(public_ontology, model)

    # filter the columns in bionty df based on fields
    bionty_df = filter_bionty_df_columns(model=model, public_ontology=public_ontology)

    # standardize in the bionty reference
    # do not inspect synonyms if the field is not name field
    inspect_synonyms = True
    if hasattr(model, "_name_field") and field.field.name != model._name_field:  # type: ignore
        inspect_synonyms = False
    result = public_ontology.inspect(
        iterable_idx,
        field=field.field.name,  # type: ignore
        mute=True,
        inspect_synonyms=inspect_synonyms,
    )
    syn_mapper = result.synonyms_mapper

    msg_syn: str = ""
    if len(syn_mapper) > 0:
        s = "" if len(syn_mapper) == 1 else "s"
        names = list(syn_mapper.keys())
        print_values = colors.purple(_format_values(names))
        msg_syn = (
            "created"
            f" {colors.purple(f'{len(syn_mapper)} {model.__name__} record{s} from Bionty')}"
            f" matching {colors.italic('synonyms')}: {print_values}"
        )

        iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # create records for values that are found in the bionty reference
    # matching either field or synonyms
    mapped_values = iterable_idx.intersection(bionty_df[field.field.name])  # type: ignore

    multi_msg = ""
    if len(mapped_values) > 0:
        bionty_kwargs, multi_msg = _bulk_create_dicts_from_df(
            keys=mapped_values,
            column_name=field.field.name,  # type: ignore
            df=bionty_df,
        )

        # # this here is needed when the organism is required to create new records
        # if organism is None:
        #     organism = _get_organism_record(
        #         field, source.organism, values=mapped_values
        #     )

        create_kwargs = (
            {"organism": organism, "source": source}
            if organism is not None
            else {"source": source}
        )
        for bk in bionty_kwargs:
            records.append(model(**bk, **create_kwargs, _skip_validation=True))

        # number of records that matches field (not synonyms)
        validated = result.validated
        if len(validated) > 0:
            s = "" if len(validated) == 1 else "s"
            print_values = colors.purple(_format_values(validated))
            # this is the success msg for existing records in the DB
            if len(msg) > 0 and not mute:
                logger.success(msg)
            if not mute:
                logger.success(
                    "created"
                    f" {colors.purple(f'{len(validated)} {model.__name__} record{s} from Bionty')}"
                    f" matching {colors.italic(f'{field.field.name}')}: {print_values}"  # type: ignore
                )

    # make sure that synonyms logging appears after the field logging
    if len(msg_syn) > 0 and not mute:
        logger.success(msg_syn)
    # warning about multi matches
    if len(multi_msg) > 0 and not mute:
        logger.warning(multi_msg)

    # return the values that are not found in the bionty reference
    unmapped_values = iterable_idx.difference(mapped_values)
    return records, unmapped_values


def index_iterable(iterable: Iterable) -> pd.Index:
    idx = pd.Index(iterable).unique()
    # No entries are made for NAs, '', None
    # returns an ordered unique not null list
    return idx[(idx != "") & (~idx.isnull())]


def _format_values(
    names: Iterable, n: int = 20, quotes: bool = True, sep: str = "'"
) -> str:
    if isinstance(names, dict):
        items = {
            f"{key}: {value}": None
            for key, value in names.items()
            if key != "None" and value != "None"
        }
    else:
        # Use a dictionary instead of a list to have unique values and preserve order
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


def _require_organism(registry: type[Record]) -> bool:
    """Check if the registry has an organism field and is required.

    Returns:
        True if the registry has an organism field and is required, False otherwise.
    """
    try:
        organism_field = registry._meta.get_field("organism")
        if organism_field.null:  # organism is not required
            return False
        else:
            return True
    except FieldDoesNotExist:
        return False


def _is_simple_field_unique(field: FieldAttr) -> bool:
    """Check if the field is an id field."""
    # id field is a unique field that's not a relation
    field = field.field
    if field.unique and not field.is_relation:
        return True
    return False


def _get_organism_record(  # type: ignore
    field: FieldAttr,
    organism: str | Record | None = None,
    values: Iterable = [],
) -> Record | None:
    """Get organism record.

    Args:
        field: the field to get the organism record for
        organism: the organism to get the record for

    Returns:
        The organism record if:
            The organism is required for the registry
            The field is not unique or the organism is not None
    """
    registry = field.field.model
    field_str = field.field.name
    check = not _is_simple_field_unique(field=field) or organism is not None

    if field_str == "ensembl_gene_id" and len(values) > 0 and organism is None:  # type: ignore
        return _organism_from_ensembl_id(values[0], registry)  # type: ignore

    if _require_organism(registry) and check:
        from bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(
            organism=organism, registry=registry, field=field_str
        )
        if organism_record is not None:
            return organism_record.save()


def _organism_from_ensembl_id(id: str) -> Record | None:  # type: ignore
    import bionty as bt

    from .artifact import Artifact  # has to be here to avoid circular imports

    ensembl_prefixes = (
        Artifact.using("laminlabs/bionty-assets")
        .get(key="ensembl_prefixes.parquet")
        .load(is_run_input=False)
        .set_index("gene_prefix")
    )
    prefix = re.sub(r"\d+", "", id)
    if prefix in ensembl_prefixes.index:
        sname = ensembl_prefixes.loc[prefix, "scientific_name"]

        organism_record = bt.Organism.from_source(scientific_name=sname)
        if organism_record is not None:
            return organism_record.save()
