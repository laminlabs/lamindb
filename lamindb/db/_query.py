import bionty as bt
import pandas as pd
from lndb_setup import settings
from sqlalchemy import inspect
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound
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


def _return_query_results_as_df(results, schema_module):
    """Return list query results as a DataFrame."""
    if len(results) > 0:
        df = pd.DataFrame([result.dict() for result in results])
    else:
        df = pd.DataFrame(columns=schema_module.__fields__.keys())

    if "id" in df.columns:
        if "v" in df.columns:
            df = df.set_index(["id", "v"])
        else:
            df = df.set_index("id")
    return df


def _create_query_func(name: str, schema_module):
    """Autogenerate query functions for each entity table."""

    def query_func(cls, as_df=False, **kwargs):
        """Query metadata from tables."""
        return Query(name=name, schema_module=schema_module, as_df=as_df, kwargs=kwargs)

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


def _get_meta_table_results(entity, link_tables, foreign_keys_backpop, entity_kwargs):
    results = getattr(query, entity)(**entity_kwargs).all()
    results_ids = [i.id for i in results]
    module_name = entity
    while module_name not in link_tables:
        if "id" not in foreign_keys_backpop[module_name]:
            return results
        parents = foreign_keys_backpop[module_name]["id"]
        for table_name, table_ref_id in parents.items():
            results = []
            for result_id in results_ids:
                results += getattr(query, table_name)(**{table_ref_id: result_id}).all()
            if table_name not in link_tables:
                results_ids = [i.id for i in results]
        module_name = table_name
    return results


def query_dobject_from_metadata(entity, entity_kwargs):
    engine = settings.instance.db_engine()
    foreign_keys = _get_all_foreign_keys(engine)
    foreign_keys_backpop = _backpopulate_foreign_keys(foreign_keys)
    link_tables = [i for i in schema.list_entities() if i.startswith("dobject_")]
    meta_results = _get_meta_table_results(
        entity=entity,
        link_tables=link_tables,
        foreign_keys_backpop=foreign_keys_backpop,
        entity_kwargs=entity_kwargs,
    )
    return meta_results


def _featureset_from_features(entity, entity_kwargs):
    """Return featuresets by quering features."""
    schema_module = getattr(schema.bionty, entity)
    stmt = _chain_select_stmt(kwargs=entity_kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")
    if len(results) > 0:
        featureset_ids = []
        for feature in results:
            featuresets = getattr(query, f"featureset_{entity}")(
                **{f"{entity}_id": feature.id}
            ).all()
            if len(featuresets) > 0:
                featureset_ids += [
                    featureset.featureset_id for featureset in featuresets
                ]
        return list(set(featureset_ids))
    return []


def query_dobject(
    id: str = None,
    v: str = None,
    name: str = None,
    dtransform_id: str = None,
    file_suffix: str = None,
    storage_id: str = None,
    time_created=None,
    time_updated=None,
    where: dict[str, dict] = None,
    as_df: bool = False,
):
    """Query from dobject."""
    schema_module = schema.core.dobject
    kwargs = {k: v for k, v in locals().items() if k in schema_module.__fields__.keys()}
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")

    if where is not None:
        dobjects = []
        for entity, entity_kwargs in where.items():
            # TODO: this part needs refactor
            try:
                bt.lookup.feature_model.__getattribute__(entity)
                # query features
                featureset_ids = _featureset_from_features(
                    entity=entity, entity_kwargs=entity_kwargs
                )
                biometas = []
                for featureset_id in featureset_ids:
                    biometas += getattr(query, "biometa")(
                        featureset_id=featureset_id
                    ).all()
                for biometa in biometas:
                    dobjects += getattr(query, "dobject_biometa")(
                        biometa_id=biometa.id
                    ).all()
            except AttributeError:
                # query obs metadata
                # find all the link tables to dobject
                dobjects += query_dobject_from_metadata(
                    entity=entity, entity_kwargs=entity_kwargs
                )

        results = [i for i in results if i.id in [j.dobject_id for j in dobjects]]

    if len(results) > 0:
        for result in results:
            track_usage(result.id, result.v, "query")

    return FilterQueryResultList(
        name="dobject", schema_module=schema_module, results=results, as_df=as_df
    )


def query_biometa(
    id: int = None,
    biosample_id: int = None,
    readout_id: int = None,
    featureset_id: int = None,
    dobject_id: str = None,
    as_df: bool = False,
):
    """Query from biometa.

    If dobject_id is provided, will search in the dobject_biometa first.
    """
    schema_module = schema.wetlab.biometa
    kwargs = {k: v for k, v in locals().items() if k in schema_module.__fields__.keys()}
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")
    # dobject_id is given, will only return results associated with dobject_id
    if dobject_id is not None:
        biometas = getattr(query, "dobject_biometa")(dobject_id=dobject_id).all()
        if len(biometas) == 0:
            results = biometas
        else:
            biometa_ids = [i.biometa_id for i in biometas]
            results = [i for i in results if i.id in biometa_ids]

    return FilterQueryResultList(
        name="biometa", schema_module=schema_module, results=results, as_df=as_df
    )


class Query:
    def __init__(self, name, schema_module, as_df=False, kwargs=None) -> None:
        self._name = name
        self._schema_module = schema_module
        self._as_df = as_df
        self._kwargs = kwargs

    def _query(self, results_type):
        stmt = _chain_select_stmt(
            kwargs=self._kwargs, schema_module=self._schema_module
        )
        results = _query_stmt(statement=stmt, results_type=results_type)
        # track usage for dobjects
        if self._name == "dobject":
            for result in results:
                track_usage(result.id, result.v, "query")
        # return DataFrame
        if self._as_df:
            return _return_query_results_as_df(
                results=results, schema_module=self._schema_module
            )
        return results

    def all(self):
        return self._query(results_type="all")

    def one(self):
        return self._query(results_type="one")

    def first(self):
        return self._query(results_type="first")


class FilterQueryResultList:
    def __init__(
        self, name: str, schema_module, results: list, as_df: bool = False
    ) -> None:
        self._name = name
        self._schema_module = schema_module
        self._results = results
        self._as_df = as_df

    def _filter(self):
        # track usage for dobjects
        if self._name == "dobject":
            for result in self._results:
                track_usage(result.id, result.v, "query")
        # return DataFrame
        if self._as_df:
            return _return_query_results_as_df(
                results=self._results, schema_module=self._schema_module
            )
        return self._results

    def all(self):
        return self._filter()

    def one(self):
        if len(self._results) == 0:
            raise NoResultFound
        elif len(self._results) > 1:
            raise MultipleResultsFound
        else:
            results = self._filter()
            if isinstance(results, pd.DataFrame):
                return results.head(1)
            return results[0]

    def first(self):
        if len(self._results) == 0:
            raise NoResultFound
        else:
            results = self._filter()
            if isinstance(results, pd.DataFrame):
                return results.head(1)
            return results[0]


class query:
    """Query literal (semantic) data."""

    pass


for name, schema_module in alltables.items():
    func = _create_query_func(name=name, schema_module=schema_module)
    setattr(query, name, classmethod(func))

setattr(query, "dobject", query_dobject)
setattr(query, "biometa", query_biometa)
