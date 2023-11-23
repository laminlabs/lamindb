from typing import Any, Dict, Iterable, List, Tuple, Union

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from lamin_utils import colors, logger
from lnschema_core.models import Feature, Registry, ULabel
from lnschema_core.types import ListLike, StrField

from .dev._settings import settings


# The base function for `from_values`
def get_or_create_records(
    iterable: ListLike,
    field: StrField,
    *,
    from_bionty: bool = False,
    **kwargs,
) -> List[Registry]:
    """Get or create records from iterables."""
    upon_create_search_names = settings.upon_create_search_names
    settings.upon_create_search_names = False
    feature: Feature = None
    try:
        Registry = field.field.model
        iterable_idx = index_iterable(iterable)

        # returns existing records & non-existing values
        records, nonexist_values, msg = get_existing_records(
            iterable_idx=iterable_idx, field=field, kwargs=kwargs
        )

        # new records to be created based on new values
        if len(nonexist_values) > 0:
            if from_bionty:
                records_bionty, unmapped_values = create_records_from_bionty(
                    iterable_idx=nonexist_values, field=field, msg=msg, **kwargs
                )
                if len(records_bionty) > 0:
                    msg = ""
                for record in records_bionty:
                    record._from_bionty = True
                records += records_bionty
            else:
                unmapped_values = nonexist_values
            # unmapped new_ids will NOT create records
            if len(unmapped_values) > 0:
                if len(msg) > 0:
                    logger.success(msg)
                s = "" if len(unmapped_values) == 1 else "s"
                print_values = colors.yellow(_print_values(unmapped_values))
                name = Registry.__name__
                n_nonval = colors.yellow(f"{len(unmapped_values)} non-validated")
                logger.warning(
                    f"{colors.red('did not create')} {name} record{s} for "
                    f"{n_nonval} {colors.italic(f'{field.field.name}{s}')}: {print_values}"  # noqa
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
    kwargs: Dict = {},
):
    model = field.field.model
    condition: Dict = {} if len(kwargs) == 0 else kwargs.copy()
    # existing records matching is agnostic to the bionty source
    if "bionty_source" in condition:
        condition.pop("bionty_source")

    if _has_organism_field(model):
        from lnschema_bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(
            organism=kwargs.get("organism"), orm=model
        )
        if organism_record is not None:
            kwargs.update({"organism": organism_record})
            condition.update({"organism": organism_record})

    # standardize based on the DB reference
    # log synonyms mapped terms
    result = model.inspect(
        iterable_idx, field=field, organism=kwargs.get("organism"), mute=True
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
        if len(msg) > 0:
            logger.success(msg)
        logger.success(syn_msg)
        msg = ""

    existing_values = iterable_idx.intersection(
        query_set.values_list(field.field.name, flat=True)
    )
    nonexist_values = iterable_idx.difference(existing_values)

    return records, nonexist_values, msg


def create_records_from_bionty(
    iterable_idx: pd.Index,
    field: StrField,
    msg: str = "",
    **kwargs,
):
    model = field.field.model
    records: List = []
    # populate additional fields from bionty
    from lnschema_bionty._bionty import get_bionty_source_record

    # create the corresponding bionty object from model
    try:
        bionty_object = model.bionty(
            organism=kwargs.get("organism"), bionty_source=kwargs.get("bionty_source")
        )
    except Exception:
        # for custom records that are not created from bionty sources
        return records, iterable_idx
    # add bionty_source record to the kwargs
    kwargs.update({"bionty_source": get_bionty_source_record(bionty_object)})

    # filter the columns in bionty df based on fields
    bionty_df = _filter_bionty_df_columns(model=model, bionty_object=bionty_object)

    # standardize in the bionty reference
    result = bionty_object.inspect(iterable_idx, field=field.field.name, mute=True)
    syn_mapper = result.synonyms_mapper

    msg_syn: str = ""
    if len(syn_mapper) > 0:
        s = "" if len(syn_mapper) == 1 else "s"
        names = list(syn_mapper.keys())
        print_values = colors.purple(_print_values(names))
        msg_syn = (
            "created"
            f" {colors.purple(f'{len(syn_mapper)} {model.__name__} record{s} from Bionty')}"  # noqa
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
        for bk in bionty_kwargs:
            records.append(model(**bk, **kwargs))

        # number of records that matches field (not synonyms)
        validated = result.validated
        if len(validated) > 0:
            s = "" if len(validated) == 1 else "s"
            print_values = colors.purple(_print_values(validated))
            # this is the success msg for existing records in the DB
            if len(msg) > 0:
                logger.success(msg)
            logger.success(
                (
                    "created"
                    f" {colors.purple(f'{len(validated)} {model.__name__} record{s} from Bionty')}"  # noqa
                    f" matching {colors.italic(f'{field.field.name}')}: {print_values}"  # noqa
                )
            )

    # make sure that synonyms logging appears after the field logging
    if len(msg_syn) > 0:
        logger.success(msg_syn)
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


def _print_values(names: List, n: int = 20) -> str:
    print_values = ", ".join([f"'{name}'" for name in names[:n]])
    if len(names) > n:
        print_values += ", ..."
    return print_values


def _filter_bionty_df_columns(model: Registry, bionty_object: Any) -> pd.DataFrame:
    bionty_df = pd.DataFrame()
    if bionty_object is not None:
        model_field_names = {i.name for i in model._meta.fields}
        # parents needs to be added here as relationships aren't in fields
        model_field_names.add("parents")
        bionty_df = bionty_object.df().reset_index()
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
    keys: Union[set, List], column_name: str, df: pd.DataFrame
) -> Tuple[Dict, str]:
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
