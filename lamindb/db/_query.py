import bionty as bt
import pandas as pd
from lndb_setup import settings
from sqlalchemy import inspect
from sqlmodel import Session, select

from .. import schema
from ..dev import track_usage
from ..schema._schema import alltables


def _query_stmt(statement, results_type="all"):
    with Session(settings.instance.db_engine()) as session:
        results = session.exec(statement).__getattribute__(results_type)()
    return results


def _chain_select_stmt(kwargs: dict, schema_module):
    stmt = select(schema_module)
    for field, value in kwargs.items():
        if (field in {"cls", "self"}) or (value is None):
            continue
        else:
            stmt = stmt.where(getattr(schema_module, field) == value)
    return stmt


def _return_query_results_as_df(statement):
    df = pd.read_sql_query(statement, settings.instance.db_engine(future=False))
    if "id" in df.columns:
        if "v" in df.columns:
            df = df.set_index(["id", "v"])
        else:
            df = df.set_index("id")
    return df


def _create_query_func(name: str, schema_module):
    def query_func(cls, return_df=False, **kwargs):
        """Query metadata from tables."""
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        if return_df:
            results = _return_query_results_as_df(statement=stmt)
        else:
            results = _query_stmt(statement=stmt, results_type="all")

        # track usage for dobjects
        if name == "dobject":
            for result in results:
                track_usage(result.id, result.v, "query")
        return results

    query_func.__name__ = name
    return query_func


def _get_all_foreign_keys(engine):
    """Result {'biosample': {'tissue_id': ('tissue', 'id')}}."""
    inspector = inspect(engine)

    def _get_foreign_keys(table_name, inspector):
        return {
            column["constrained_columns"][0]: (
                column["referred_table"],
                column["referred_columns"][0],
            )
            for column in inspector.get_foreign_keys(table_name)
        }

    foreign_keys = {}
    for table_name in inspector.get_table_names():
        foreign_keys_table = _get_foreign_keys(table_name, inspector)
        if len(foreign_keys_table) > 0:
            foreign_keys[table_name] = foreign_keys_table

    return foreign_keys


def _backpopulate_foreign_keys(foreign_keys):
    """Result {'tissue': {'id': {'biosample': 'tissue_id'}}}."""
    foreign_keys_backpop = {}

    for module_name, keys in foreign_keys.items():
        for key, (module, ref_key) in keys.items():
            if foreign_keys_backpop.get(module) is None:
                foreign_keys_backpop[module] = {}
            if foreign_keys_backpop[module].get(ref_key) is None:
                foreign_keys_backpop[module][ref_key] = {}
            foreign_keys_backpop[module][ref_key][module_name] = key

    return foreign_keys_backpop


def _get_meta_table_results(
    entity_name, link_tables, foreign_keys_backpop, **entity_kwargs
):
    results = getattr(query, entity_name)(**entity_kwargs)
    results_ids = [i.id for i in results]
    module_name = entity_name
    while module_name not in link_tables:
        if "id" not in foreign_keys_backpop[module_name]:
            return results
        parents = foreign_keys_backpop[module_name]["id"]
        for table_name, table_ref_id in parents.items():
            results = []
            for result_id in results_ids:
                results += getattr(query, table_name)(**{table_ref_id: result_id})
            if table_name not in link_tables:
                results_ids = [i.id for i in results]
        module_name = table_name
    return results


def query_dobject_from_metadata(entity_name, **entity_kwargs):
    engine = settings.instance.db_engine()
    foreign_keys = _get_all_foreign_keys(engine)
    foreign_keys_backpop = _backpopulate_foreign_keys(foreign_keys)
    link_tables = [i for i in schema.list_entities() if i.startswith("dobject_")]
    meta_results = _get_meta_table_results(
        entity_name=entity_name,
        link_tables=link_tables,
        foreign_keys_backpop=foreign_keys_backpop**entity_kwargs,
    )
    dobject_ids = set([dobject.dobject_id for dobject in meta_results])
    if len(dobject_ids) > 0:
        dobjects = []
        for dobject_id in dobject_ids:
            dobjects += getattr(query, "dobject")(id=dobject_id)
        return dobjects
    return []


class query:
    """Query literal (semantic) data."""

    @classmethod
    def table_as_df(cls, entity_name) -> pd.DataFrame:
        """Load metadata table as DataFrame."""
        engine = settings.instance.db_engine()
        with engine.connect() as conn:
            df = pd.read_sql_table(entity_name, conn)
            if "id" in df.columns:
                if "v" in df.columns:
                    df = df.set_index(["id", "v"])
                else:
                    df = df.set_index("id")
        return df


def _featureset_from_features(entity_name, **entity_kwargs):
    """Return featuresets by quering features."""
    schema_module = getattr(schema.bionty, entity_name)
    stmt = _chain_select_stmt(kwargs=entity_kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")
    if len(results) > 0:
        featureset_ids = []
        for feature in results:
            featuresets = getattr(query, f"featureset_{entity_name}")(
                **{f"{entity_name}_id": feature.id}
            )
            if len(featuresets) > 0:
                featureset_ids += [
                    featureset.featureset_id for featureset in featuresets
                ]
        return list(set(featureset_ids))
    return []


def dobject(
    id: str = None,
    v: str = None,
    name: str = None,
    dtransform_id: str = None,
    file_suffix: str = None,
    storage_id: str = None,
    time_created=None,
    time_updated=None,
    entity_name: str = None,
    **entity_kwargs,
):
    """Query from dobject."""
    schema_module = schema.core.dobject
    kwargs = {
        k: v for k, v in locals().items() if k in schema.core.dobject.__fields__.keys()
    }
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")

    if entity_name is not None:
        if getattr(bt.lookup.feature_model, entity_name) is not None:
            # query features
            featureset_ids = _featureset_from_features(
                entity_name=entity_name, **entity_kwargs
            )
            biometas = []
            for featureset_id in featureset_ids:
                biometas += getattr(query, "biometa")(featureset_id=featureset_id)
            dobjects = []
            for biometa in biometas:
                dobjects += getattr(query, "dobject_biometa")(biometa_id=biometa.id)
        else:
            # query obs metadata
            # find all the link tables to dobject
            results = query_dobject_from_metadata(
                entity_name=entity_name, **entity_kwargs
            )

        results = [i for i in results if i.id in [j.dobject_id for j in dobjects]]

    if len(results) > 0:
        for result in results:
            track_usage(result.id, result.v, "query")

    return results


def biometa(
    id: int = None,
    biosample_id: int = None,
    readout_id: int = None,
    featureset_id: int = None,
    dobject_id: str = None,
):
    """Query from biometa.

    If dobject_id is provided, will search in the dobject_biometa first.
    """
    kwargs = locals()
    del kwargs["dobject_id"]
    schema_module = schema.wetlab.biometa
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")
    # dobject_id is given, will only return results associated with dobject_id
    if dobject_id is not None:
        biometas = getattr(query, "dobject_biometa")(dobject_id=dobject_id)
        if len(biometas) == 0:
            return biometas
        else:
            biometa_ids = [i.biometa_id for i in biometas]
            return [i for i in results if i.id in biometa_ids]
    else:
        return


for name, schema_module in alltables.items():
    func = _create_query_func(name=name, schema_module=schema_module)
    setattr(query, name, classmethod(func))

setattr(query, "dobject", dobject)
setattr(query, "biometa", biometa)
