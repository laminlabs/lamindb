from __future__ import annotations

from typing import TYPE_CHECKING, Iterable

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from lamin_utils import colors, logger
from lnschema_core.models import Feature, Record, ULabel

from .core._settings import settings

if TYPE_CHECKING:
    from lnschema_core.types import ListLike, StrField


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
) -> list[Record]:
    """Get or create records from iterables."""
    Record = field.field.model
    if create:
        return [Record(**{field.field.name: value}) for value in iterable]
    creation_search_names = settings.creation.search_names
    feature: Feature = None
    organism = _get_organism_record(field, organism)
    kwargs: dict = {}
    if organism is not None:
        kwargs["organism"] = organism
    if source is not None:
        kwargs["source"] = source
    settings.creation.search_names = False
    try:
        iterable_idx = index_iterable(iterable)

        # returns existing records & non-existing values
        records, nonexist_values, msg = get_existing_records(
            iterable_idx=iterable_idx, field=field, mute=mute, **kwargs
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
            if not source_record and hasattr(Record, "public"):
                from bionty._bionty import get_source_record

                source_record = get_source_record(Record.public(organism=organism))
            if source_record:
                from bionty.core._add_ontology import check_source_in_db

                check_source_in_db(
                    registry=Record,
                    source=source_record,
                    update=True,
                )

                from_source = not source_record.in_db
            elif hasattr(Record, "source_id"):
                from_source = True
            else:
                from_source = False

            if from_source:
                records_bionty, unmapped_values = create_records_from_source(
                    iterable_idx=nonexist_values,
                    field=field,
                    msg=msg,
                    mute=mute,
                    **kwargs,
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
                print_values = colors.yellow(_print_values(unmapped_values))
                name = Record.__name__
                n_nonval = colors.yellow(f"{len(unmapped_values)} non-validated")
                if not mute:
                    logger.warning(
                        f"{colors.red('did not create')} {name} record{s} for "
                        f"{n_nonval} {colors.italic(f'{field.field.name}{s}')}: {print_values}"
                    )
        if Record.__module__.startswith("bionty.") or Record == ULabel:
            if isinstance(iterable, pd.Series):
                feature = iterable.name
            feature_name = None
            if isinstance(feature, str):
                feature_name = feature
            if feature_name is not None:
                if feature_name is not None:
                    for record in records:
                        record._feature = feature_name
                logger.debug(f"added default feature '{feature_name}'")
        return records
    finally:
        settings.creation.search_names = creation_search_names


def get_existing_records(
    iterable_idx: pd.Index,
    field: StrField,
    mute: bool = False,
    **kwargs,
):
    model = field.field.model
    condition: dict = {} if len(kwargs) == 0 else kwargs.copy()
    # existing records matching is agnostic to the bionty source
    if "source" in condition:
        condition.pop("source")

    # standardize based on the DB reference
    # log synonyms mapped terms
    result = model.inspect(
        iterable_idx,
        field=field,
        organism=kwargs.get("organism"),
        source=kwargs.get("source"),
        mute=True,
    )
    syn_mapper = result.synonyms_mapper

    syn_msg = ""
    if len(syn_mapper) > 0:
        s = "" if len(syn_mapper) == 1 else "s"
        names = list(syn_mapper.keys())
        print_values = colors.green(_print_values(names))
        syn_msg = (
            "loaded"
            f" {colors.green(f'{len(syn_mapper)} {model.__name__} record{s}')}"
            f" matching {colors.italic('synonyms')}: {print_values}"
        )
        iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # get all existing records in the db
    # if necessary, create records for the values in kwargs
    # k:v -> k:v_record
    # kwargs is used to deal with organism
    condition.update({f"{field.field.name}__in": iterable_idx.values})

    query_set = model.filter(**condition)
    records = query_set.list()

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
    validated = result.validated
    msg = ""
    if len(validated) > 0:
        s = "" if len(validated) == 1 else "s"
        print_values = colors.green(_print_values(validated))
        msg = (
            "loaded"
            f" {colors.green(f'{len(validated)} {model.__name__} record{s}')}"
            f" matching {colors.italic(f'{field.field.name}')}: {print_values}"
        )

    # no logging if all values are validated
    # logs if there are synonyms
    if len(syn_msg) > 0:
        if len(msg) > 0 and not mute:
            logger.success(msg)
        if not mute:
            logger.success(syn_msg)
        msg = ""

    existing_values = iterable_idx.intersection(
        query_set.values_list(field.field.name, flat=True)
    )
    nonexist_values = iterable_idx.difference(existing_values)

    return records, nonexist_values, msg


def create_records_from_source(
    iterable_idx: pd.Index,
    field: StrField,
    msg: str = "",
    mute: bool = False,
    **kwargs,
):
    model = field.field.model
    records: list = []
    # populate additional fields from bionty
    from bionty._bionty import get_source_record
    from bionty.core._bionty import filter_bionty_df_columns

    # create the corresponding bionty object from model
    try:
        # TODO: more generic
        organism = kwargs.get("organism")
        if field.field.name == "ensembl_gene_id":
            if iterable_idx[0].startswith("ENSG"):
                organism = "human"
            elif iterable_idx[0].startswith("ENSMUSG"):
                organism = "mouse"
        public_ontology = model.public(organism=organism, source=kwargs.get("source"))
    except Exception:
        # for custom records that are not created from public sources
        return records, iterable_idx
    # add source record to the kwargs
    source_record = get_source_record(public_ontology)
    kwargs.update({"source": source_record})

    # filter the columns in bionty df based on fields
    bionty_df = filter_bionty_df_columns(model=model, public_ontology=public_ontology)

    # standardize in the bionty reference
    result = public_ontology.inspect(iterable_idx, field=field.field.name, mute=True)
    syn_mapper = result.synonyms_mapper

    msg_syn: str = ""
    if len(syn_mapper) > 0:
        s = "" if len(syn_mapper) == 1 else "s"
        names = list(syn_mapper.keys())
        print_values = colors.purple(_print_values(names))
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
        organism_kwargs = {}
        if "organism" not in kwargs:
            organism_record = _get_organism_record(
                field, public_ontology.organism, force=True
            )
            if organism_record is not None:
                organism_kwargs["organism"] = organism_record
        for bk in bionty_kwargs:
            records.append(model(**bk, **kwargs, **organism_kwargs))

        # number of records that matches field (not synonyms)
        validated = result.validated
        if len(validated) > 0:
            s = "" if len(validated) == 1 else "s"
            print_values = colors.purple(_print_values(validated))
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


def _print_values(names: Iterable, n: int = 20, quotes: bool = True) -> str:
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
        unique_items = [f"'{item}'" for item in unique_items]

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
            print_values = _print_values(dup)
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
    registry = field.field.model
    check = True
    if not force and hasattr(registry, "_ontology_id_field"):
        check = field.field.name != registry._ontology_id_field
        # e.g. bionty.CellMarker has "name" as _ontology_id_field
        if not registry._ontology_id_field.endswith("id"):
            check = True

    if _has_organism_field(registry) and check:
        from bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(organism=organism, orm=registry)
        if organism_record is not None:
            return organism_record
