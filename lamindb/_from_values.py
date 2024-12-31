from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from lamin_utils import colors, logger

from lamindb._query_set import RecordList
from lamindb.models import Record

from .core._settings import settings

if TYPE_CHECKING:
    from collections.abc import Iterable

    from lamindb.types import ListLike, StrField


# The base function for `from_values`
def get_or_create_records(
    iterable: ListLike,
    field: StrField,
    *,
    create: bool = False,
    from_source: bool = False,
    organism: Record | str | None = None,
    source: Record | None = None,
    mute: bool = False,
) -> RecordList:
    """Get or create records from iterables."""
    registry = field.field.model
    if create:
        return RecordList([registry(**{field.field.name: value}) for value in iterable])
    creation_search_names = settings.creation.search_names
    organism = _get_organism_record(field, organism)
    settings.creation.search_names = False
    try:
        iterable_idx = index_iterable(iterable)

        # returns existing records & non-existing values
        records, nonexist_values, msg = get_existing_records(
            iterable_idx=iterable_idx,
            field=field,
            organism=organism,
            mute=mute,
        )

        # new records to be created based on new values
        if len(nonexist_values) > 0:
            source_record = None
            if from_source:
                if isinstance(source, Record):
                    source_record = source
                elif (
                    len(records) > 0
                    and hasattr(records[0], "source_id")
                    and records[0].source_id
                ):
                    source_record = records[0].source
            if not source_record and hasattr(registry, "public"):
                if organism is None:
                    organism = _ensembl_prefix(nonexist_values[0], field, organism)
                    organism = _get_organism_record(field, organism, force=True)

            if source_record:
                from bionty.core._add_ontology import check_source_in_db

                check_source_in_db(registry=registry, source=source_record)

                from_source = not source_record.in_db
            elif hasattr(registry, "source_id"):
                from_source = True
            else:
                from_source = False

            if from_source:
                records_bionty, unmapped_values = create_records_from_source(
                    iterable_idx=nonexist_values,
                    field=field,
                    organism=organism,
                    source=source_record,
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
                if len(msg) > 0 and not mute:
                    logger.success(msg)
                s = "" if len(unmapped_values) == 1 else "s"
                print_values = colors.yellow(_format_values(unmapped_values))
                name = registry.__name__
                n_nonval = colors.yellow(f"{len(unmapped_values)} non-validated")
                if not mute:
                    logger.warning(
                        f"{colors.red('did not create')} {name} record{s} for "
                        f"{n_nonval} {colors.italic(f'{field.field.name}{s}')}: {print_values}"
                    )
        # if registry.__get_schema_name__() == "bionty" or registry == ULabel:
        #     if isinstance(iterable, pd.Series):
        #         feature = iterable.name
        #     feature_name = None
        #     if isinstance(feature, str):
        #         feature_name = feature
        #     if feature_name is not None:
        #         if feature_name is not None:
        #             for record in records:
        #                 record._feature = feature_name
        #         logger.debug(f"added default feature '{feature_name}'")
        return RecordList(records)
    finally:
        settings.creation.search_names = creation_search_names


def get_existing_records(
    iterable_idx: pd.Index,
    field: StrField,
    organism: Record | None = None,
    mute: bool = False,
):
    # NOTE: existing records matching is agnostic to the source
    model = field.field.model
    if organism is None and field.field.name == "ensembl_gene_id":
        if len(iterable_idx) > 0:
            organism = _ensembl_prefix(iterable_idx[0], field, organism)
            organism = _get_organism_record(field, organism, force=True)

    # standardize based on the DB reference
    # log synonyms mapped terms
    syn_mapper = model.standardize(
        iterable_idx,
        field=field,
        organism=organism,
        mute=True,
        public_aware=False,
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
    # if necessary, create records for the values in kwargs
    # k:v -> k:v_record
    query = {f"{field.field.name}__in": iterable_idx.values}
    if organism is not None:
        query["organism"] = organism
    records = model.filter(**query).list()

    if len(validated) == len(iterable_idx):
        return records, [], msg
    else:
        nonval_values = iterable_idx.difference(validated)
        return records, nonval_values, msg


def create_records_from_source(
    iterable_idx: pd.Index,
    field: StrField,
    organism: Record | None = None,
    source: Record | None = None,
    msg: str = "",
    mute: bool = False,
):
    model = field.field.model
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
    result = public_ontology.inspect(iterable_idx, field=field.field.name, mute=True)
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
    mapped_values = iterable_idx.intersection(bionty_df[field.field.name])

    multi_msg = ""
    if len(mapped_values) > 0:
        bionty_kwargs, multi_msg = _bulk_create_dicts_from_df(
            keys=mapped_values, column_name=field.field.name, df=bionty_df
        )

        if hasattr(model, "organism_id") and organism is None:
            organism = _get_organism_record(field, source.organism, force=True)

        create_kwargs = (
            {"organism": organism, "source": source}
            if organism is not None
            else {"source": source}
        )
        for bk in bionty_kwargs:
            records.append(model(**bk, **create_kwargs))

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
                    f" matching {colors.italic(f'{field.field.name}')}: {print_values}"
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


def _has_organism_field(registry: type[Record]) -> bool:
    try:
        registry._meta.get_field("organism")
        return True
    except FieldDoesNotExist:
        return False


def _get_organism_record(
    field: StrField, organism: str | Record, force: bool = False
) -> Record:
    """Get organism record.

    Args:
        field: the field to get the organism record for
        organism: the organism to get the record for
        force: whether to force fetching the organism record
    """
    registry = field.field.model
    check = True
    if not force and hasattr(registry, "_ontology_id_field"):
        check = field.field.name != registry._ontology_id_field
        # e.g. bionty.CellMarker has "name" as _ontology_id_field
        if not registry._ontology_id_field.endswith("id"):
            check = True

    if _has_organism_field(registry) and check:
        from bionty._bionty import create_or_get_organism_record

        if field and not isinstance(field, str):
            field = field.field.name

        organism_record = create_or_get_organism_record(
            organism=organism, registry=registry, field=field
        )
        if organism_record is not None:
            return organism_record


def _ensembl_prefix(id: str, field: StrField, organism: Record | None) -> str | None:
    if field.field.name == "ensembl_gene_id" and organism is None:
        if id.startswith("ENSG"):
            organism = "human"
        elif id.startswith("ENSMUSG"):
            organism = "mouse"

    return organism
