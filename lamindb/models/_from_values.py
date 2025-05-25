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
    mute: bool = False,
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
    )

    # new records to be created based on new values
    if len(nonexist_values) > 0:
        if registry.__base__.__name__ == "BioRecord":
            from bionty._organism import is_organism_required

            # if can and needed, get organism record from the existing records
            if (
                organism_record is None
                and len(records) > 0
                and is_organism_required(registry)
            ):
                organism_record = records[0].organism
            records_public, unmapped_values = create_records_from_source(
                iterable_idx=nonexist_values,
                field=field,
                organism=organism_record,
                source=source,
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
    mute: bool = False,
) -> tuple[list, pd.Index, str]:
    """Get existing records from the database."""
    # NOTE: existing records matching is agnostic to the source
    model = field.field.model  # type: ignore

    # log synonyms mapped terms
    syn_mapper = model.standardize(
        iterable_idx,
        field=field,
        organism=organism,
        mute=True,
        source_aware=False,  # standardize only based on the DB reference
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
    organism: SQLRecord | None = None,
    source: SQLRecord | None = None,
    msg: str = "",
    mute: bool = False,
) -> tuple[list, pd.Index]:
    """Create records from source."""
    model = field.field.model  # type: ignore
    records: list = []
    # populate additional fields from public_df
    from bionty._source import filter_public_df_columns, get_source_record

    # get the default source
    source_record = get_source_record(model, organism, source)

    # create the corresponding PublicOntology object from model
    try:
        public_ontology = model.public(source=source_record)
    except Exception:
        # no public source
        return records, iterable_idx

    # filter the columns in public df based on fields
    public_df = filter_public_df_columns(model=model, public_ontology=public_ontology)

    if public_df.empty:
        return records, iterable_idx

    # standardize in the public reference
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

        # this here is needed when the organism is required to create new records
        if organism is None:
            organism = get_organism_record_from_field(
                field, source_record.organism, values=mapped_values
            )

        create_kwargs = (
            {"organism": organism, "source": source_record}
            if organism is not None
            else {"source": source_record}
        )
        for bk in public_kwargs:
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


def get_organism_record_from_field(  # type: ignore
    field: FieldAttr,
    organism: str | SQLRecord | None = None,
    values: ListLike = None,
    using_key: str | None = None,
) -> SQLRecord | None:
    """Get organism record.

    Args:
        field: the field to get the organism record for
        organism: the organism to get the record for
        values: the values to get the organism record for
        using_key: the db to get the organism record for

    Returns:
        The organism record if:
            The organism FK is required for the registry
            The field is not unique or the organism is not None
    """
    if values is None:
        values = []
    registry = field.field.model
    field_str = field.field.name
    # id field is a unique field that's not a relation
    is_simple_field_unique = field.field.unique and not field.field.is_relation
    check = not is_simple_field_unique or organism is not None

    if (
        registry.__get_name_with_module__() == "bionty.Gene"
        and field.field.name == "ensembl_gene_id"
        and len(values) > 0
        and organism is None
    ):  # type: ignore
        from bionty._organism import organism_from_ensembl_id

        return organism_from_ensembl_id(values[0], using_key)  # type: ignore

    if registry.__base__.__name__ == "BioRecord" and check:
        from bionty._organism import create_or_get_organism_record

        organism_record = create_or_get_organism_record(
            organism=organism, registry=registry, field=field_str
        )
        return organism_record
