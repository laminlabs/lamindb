from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_utils import colors, logger
from lnschema_core.models import ORM, Feature
from lnschema_core.types import ListLike

from .dev._settings import settings


# The base function for `from_values`
def get_or_create_records(
    iterable: ListLike,
    field: Field,
    *,
    from_bionty: bool = False,
    **kwargs,
) -> List:
    """Get or create records from iterables."""
    upon_create_search_names = settings.upon_create_search_names
    settings.upon_create_search_names = False
    feature: Feature = None
    if "feature" in kwargs:
        feature = kwargs.pop("feature")
        kwargs["feature_id"] = feature.id
    types: Optional[Dict] = None
    if "types" in kwargs:
        types = kwargs.pop("types")
    try:
        field_name = field.field.name
        ORM = field.field.model
        iterable_idx = index_iterable(iterable)

        if isinstance(ORM, Feature):
            if types is None:
                raise ValueError("Please pass types as {} or use FeatureSet.from_df()")

        # returns existing records & non-existing values
        records, nonexist_values = get_existing_records(
            iterable_idx=iterable_idx, field=field, kwargs=kwargs
        )

        # new records to be created based on new values
        if len(nonexist_values) > 0:
            if from_bionty:
                records_bionty, unmapped_values = create_records_from_bionty(
                    iterable_idx=nonexist_values, field=field, **kwargs
                )
                records += records_bionty
            else:
                unmapped_values = nonexist_values
            # unmapped new_ids will only create records with field and kwargs
            if len(unmapped_values) > 0:
                for value in unmapped_values:
                    params = {field_name: value}
                    if types is not None:
                        params["type"] = str(types[value])
                    records.append(ORM(**params, **kwargs))
                s = "" if len(unmapped_values) == 1 else "s"
                print_unmapped_values = ", ".join(unmapped_values[:7])
                if len(unmapped_values) > 10:
                    print_unmapped_values += ", ..."
                additional_info = " "
                if feature is not None:
                    additional_info = f" Feature {feature.name} and "
                logger.warning(
                    f"Created {colors.yellow(f'{len(unmapped_values)} {ORM.__name__} record{s}')} for{additional_info}"  # noqa
                    f"{colors.yellow(f'{field_name}{s}')}: {print_unmapped_values}"  # noqa
                )
        return records
    finally:
        settings.upon_create_search_names = upon_create_search_names


def get_existing_records(iterable_idx: pd.Index, field: Field, kwargs: Dict = {}):
    field_name = field.field.name
    model = field.field.model
    condition: Dict = {}

    if _has_species_field(model):
        from lnschema_bionty._bionty import create_or_get_species_record

        species_record = create_or_get_species_record(
            species=kwargs.get("species"), orm=model
        )
        if species_record is not None:
            kwargs.update({"species": species_record})
            condition.update({"species__name": species_record.name})

    # map synonyms based on the DB reference
    syn_mapper = model.map_synonyms(
        iterable_idx, species=kwargs.get("species"), return_mapper=True
    )

    syn_msg = ""
    if len(syn_mapper) > 0:
        s = "" if len(syn_mapper) == 1 else "s"
        syn_msg = (
            "Loaded"
            f" {colors.green(f'{len(syn_mapper)} {model.__name__} record{s}')} that"  # noqa
            f" matched {colors.green('synonyms')}"
        )
        iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # get all existing records in the db
    # if necessary, create records for the values in kwargs
    # k:v -> k:v_record
    # kwargs is used to deal with species
    condition.update({f"{field_name}__in": iterable_idx.values})

    from ._select import select

    stmt = select(model, **condition)

    records = stmt.list()  # existing records
    n_name = len(records) - len(syn_mapper)
    if n_name > 0:
        s = "" if n_name == 1 else "s"
        logger.info(
            "Loaded"
            f" {colors.green(f'{n_name} {model.__name__} record{s}')} that"
            f" matched {colors.green(f'{field_name}')}"
        )
    # make sure that synonyms logging appears after the field logging
    if len(syn_msg) > 0:
        logger.info(syn_msg)

    existing_values = iterable_idx.intersection(stmt.values_list(field_name, flat=True))
    nonexist_values = iterable_idx.difference(existing_values)

    return records, nonexist_values


def create_records_from_bionty(
    iterable_idx: pd.Index,
    field: Field,
    **kwargs,
):
    model = field.field.model
    field_name = field.field.name
    records: List = []
    # populate additional fields from bionty
    from lnschema_bionty._bionty import get_bionty_source_record

    # create the corresponding bionty object from model
    bionty_object = model.bionty(species=kwargs.get("species"))
    # add bionty_source record to the kwargs
    kwargs.update({"bionty_source": get_bionty_source_record(bionty_object)})

    # filter the columns in bionty df based on fields
    bionty_df = _filter_bionty_df_columns(model=model, bionty_object=bionty_object)

    # map synonyms in the bionty reference
    try:
        syn_mapper = bionty_object.map_synonyms(iterable_idx, return_mapper=True)
    except KeyError:
        # no synonyms column
        syn_mapper = {}
    msg_syn: str = ""
    if len(syn_mapper) > 0:
        s = "" if len(syn_mapper) == 1 else "s"
        msg_syn = (
            "Loaded"
            f" {colors.purple(f'{len(syn_mapper)} {model.__name__} record{s} from Bionty')} that"  # noqa
            f" matched {colors.purple('synonyms')}"
        )

        iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # create records for values that are found in the bionty reference
    mapped_values = iterable_idx.intersection(bionty_df[field_name])

    if len(mapped_values) > 0:
        bionty_kwargs, multi_msg = _bulk_create_dicts_from_df(
            keys=mapped_values, column_name=field_name, df=bionty_df
        )
        for bk in bionty_kwargs:
            records.append(model(**bk, **kwargs))

        # logging of BiontySource linking
        source_msg = (
            ""
            if kwargs.get("bionty_source") is None
            else f" (bionty_source_id={kwargs.get('bionty_source').id})"  # type:ignore # noqa
        )

        # number of records that matches field (not synonyms)
        n_name = len(records) - len(syn_mapper)
        if n_name > 0:
            s = "" if n_name == 1 else "s"
            msg = (
                "Loaded"
                f" {colors.purple(f'{n_name} {model.__name__} record{s} from Bionty')} that"  # noqa
                f" matched {colors.purple(f'{field_name}')}"
            )
            logger.info(msg + source_msg)
        # make sure that synonyms logging appears after the field logging
        if len(msg_syn) > 0:
            logger.info(msg_syn + source_msg)
        # warning about multi matches
        if len(multi_msg) > 0:
            logger.warning(multi_msg)

    # return the values that are not found in the bionty reference
    unmapped_values = iterable_idx.difference(mapped_values)
    return records, unmapped_values


def index_iterable(iterable: Iterable) -> pd.Index:
    idx = pd.Index(iterable).unique()
    # No entries are made for NAs, '', None
    # returns an ordered unique not null list
    return idx[(idx != "") & (~idx.isnull())]


def _filter_bionty_df_columns(model: ORM, bionty_object: Any) -> pd.DataFrame:
    bionty_df = pd.DataFrame()
    if bionty_object is not None:
        model_field_names = {i.name for i in model._meta.fields}
        # parents needs to be added here as relationships aren't in fields
        model_field_names.add("parents")
        bionty_df = bionty_object.df().reset_index()
        if model.__name__ == "Gene":
            # groupby ensembl_gene_id and concat ncbi_gene_ids
            bionty_df.drop(
                columns=["hgnc_id", "mgi_id", "index"], errors="ignore", inplace=True
            )
            bionty_df.drop_duplicates(["ensembl_gene_id", "ncbi_gene_id"], inplace=True)
            bionty_df["ncbi_gene_id"] = bionty_df["ncbi_gene_id"].fillna("")
            bionty_df = (
                bionty_df.groupby("ensembl_gene_id")
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
    keys: Union[set, List], column_name: str, df: pd.DataFrame
) -> Tuple[Dict, str]:
    """Get fields from a DataFrame for many rows."""
    multi_msg = ""
    if df.index.name != column_name:
        df = df.set_index(column_name).loc[list(keys)]
    if not df.index.is_unique:
        # return all records for multi-matches with a warning
        dup = df.index[df.index.duplicated()].unique().tolist()
        multi_msg = f"Multiple matches found in Bionty for: {dup}"

    return df.reset_index().to_dict(orient="records"), multi_msg


def _has_species_field(orm: ORM) -> bool:
    try:
        orm._meta.get_field("species")
        return True
    except FieldDoesNotExist:
        return False
