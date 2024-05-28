from __future__ import annotations

from typing import TYPE_CHECKING, Any, Iterable

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from lamin_utils import colors, logger
from lnschema_core.models import Feature, Registry, ULabel

from .core._settings import settings

if TYPE_CHECKING:
    from lnschema_core.types import ListLike, StrField


# The base function for `from_values`
def get_or_create_records(
    iterable: ListLike,
    field: StrField,
    *,
    create: bool = False,
    from_public: bool = False,
    organism: Registry | str | None = None,
    public_source: Registry | None = None,
    mute: bool = False,
) -> list[Registry]:
    """Get or create records from iterables."""
    Registry = field.field.model
    if create:
        return [Registry(**{field.field.name: value}) for value in iterable]
    upon_create_search_names = settings.upon_create_search_names
    feature: Feature = None
    organism = _get_organism_record(field, organism)
    kwargs: dict = {}
    if organism is not None:
        kwargs["organism"] = organism
    if public_source is not None:
        kwargs["public_source"] = public_source
    settings.upon_create_search_names = False
    try:
        iterable_idx = index_iterable(iterable)

        # returns existing records & non-existing values
        records, nonexist_values, msg = get_existing_records(
            iterable_idx=iterable_idx, field=field, mute=mute, **kwargs
        )

        # new records to be created based on new values
        if len(nonexist_values) > 0:
            if from_public:
                records_bionty, unmapped_values = create_records_from_public(
                    iterable_idx=nonexist_values,
                    field=field,
                    msg=msg,
                    mute=mute,
                    **kwargs,
                )
                if len(records_bionty) > 0:
                    msg = ""
                for record in records_bionty:
                    record._from_public = True
                records += records_bionty
            else:
                unmapped_values = nonexist_values
            # unmapped new_ids will NOT create records
            if len(unmapped_values) > 0:
                if len(msg) > 0 and not mute:
                    logger.success(msg)
                s = "" if len(unmapped_values) == 1 else "s"
                print_values = colors.yellow(_print_values(unmapped_values))
                name = Registry.__name__
                n_nonval = colors.yellow(f"{len(unmapped_values)} non-validated")
                if not mute:
                    logger.warning(
                        f"{colors.red('did not create')} {name} record{s} for "
                        f"{n_nonval} {colors.italic(f'{field.field.name}{s}')}: {print_values}"
                    )
        if Registry.__module__.startswith("lnschema_bionty.") or Registry == ULabel:
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
        settings.upon_create_search_names = upon_create_search_names


def get_existing_records(
    iterable_idx: pd.Index,
    field: StrField,
    mute: bool = False,
    **kwargs,
):
    model = field.field.model
    condition: dict = {} if len(kwargs) == 0 else kwargs.copy()
    # existing records matching is agnostic to the bionty source
    if "public_source" in condition:
        condition.pop("public_source")

    # standardize based on the DB reference
    # log synonyms mapped terms
    result = model.inspect(
        iterable_idx,
        field=field,
        organism=kwargs.get("organism"),
        public_source=kwargs.get("public_source"),
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


def create_records_from_public(
    iterable_idx: pd.Index,
    field: StrField,
    msg: str = "",
    mute: bool = False,
    **kwargs,
):
    model = field.field.model
    records: list = []
    # populate additional fields from bionty
    from lnschema_bionty._bionty import get_public_source_record

    # create the corresponding bionty object from model
    try:
        # TODO: more generic
        organism = kwargs.get("organism")
        if field.field.name == "ensembl_gene_id":
            if iterable_idx[0].startswith("ENSG"):
                organism = "human"
            elif iterable_idx[0].startswith("ENSMUSG"):
                organism = "mouse"
        public_ontology = model.public(
            organism=organism, public_source=kwargs.get("public_source")
        )
    except Exception:
        # for custom records that are not created from public sources
        return records, iterable_idx
    # add public_source record to the kwargs
    kwargs.update({"public_source": get_public_source_record(public_ontology)})

    # filter the columns in bionty df based on fields
    bionty_df = _filter_bionty_df_columns(model=model, public_ontology=public_ontology)

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
    names = (name for name in names if name != "None")
    unique_names = list(dict.fromkeys(names))[:n]
    if quotes:
        print_values = ", ".join(f"'{name}'" for name in unique_names)
    else:
        print_values = ", ".join(f"{name}" for name in unique_names)
    if len(unique_names) > n:
        print_values += ", ..."
    return print_values


def _filter_bionty_df_columns(model: Registry, public_ontology: Any) -> pd.DataFrame:
    bionty_df = pd.DataFrame()
    if public_ontology is not None:
        model_field_names = {i.name for i in model._meta.fields}
        # parents needs to be added here as relationships aren't in fields
        model_field_names.add("parents")
        bionty_df = public_ontology.df().reset_index()
        if model.__name__ == "Gene":
            # groupby ensembl_gene_id and concat ncbi_gene_ids
            groupby_id_col = (
                "ensembl_gene_id" if "ensembl_gene_id" in bionty_df else "stable_id"
            )
            bionty_df.drop(
                columns=["hgnc_id", "mgi_id", "index"], errors="ignore", inplace=True
            )
            bionty_df.drop_duplicates([groupby_id_col, "ncbi_gene_id"], inplace=True)
            bionty_df["ncbi_gene_id"] = bionty_df["ncbi_gene_id"].fillna("")
            bionty_df = (
                bionty_df.groupby(groupby_id_col)
                .agg(
                    {
                        "symbol": "first",
                        "ncbi_gene_id": "|".join,
                        "biotype": "first",
                        "description": "first",
                        "synonyms": "first",
                    }
                )
                .reset_index()
            )
            bionty_df.rename(columns={"ncbi_gene_id": "ncbi_gene_ids"}, inplace=True)
        # rename definition to description for the lnschema_bionty
        bionty_df.rename(columns={"definition": "description"}, inplace=True)
        bionty_df = bionty_df.loc[:, bionty_df.columns.isin(model_field_names)]
    return bionty_df


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


def _has_organism_field(orm: Registry) -> bool:
    try:
        orm._meta.get_field("organism")
        return True
    except FieldDoesNotExist:
        return False


def _get_organism_record(
    field: StrField, organism: str | Registry, force: bool = False
) -> Registry:
    model = field.field.model
    check = True if force else field.field.name != "ensembl_gene_id"

    if _has_organism_field(model) and check:
        from lnschema_bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(organism=organism, orm=model)
        if organism_record is not None:
            return organism_record
