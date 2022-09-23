from functools import cached_property
from typing import Dict

import bionty as bt
import pandas as pd
from lndb_setup import settings
from sqlalchemy import inspect
from sqlmodel import Session, select

from ..dev import track_usage
from ..dev.db import exception
from ..schema._table import Table


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


def _return_query_results_as_df(results, model):
    """Return list query results as a DataFrame."""
    if len(results) > 0:
        df = pd.DataFrame([result.dict() for result in results])
    else:
        df = pd.DataFrame(columns=Table.get_fields(model))

    if "id" in df.columns:
        if "v" in df.columns:
            df = df.set_index(["id", "v"])
        else:
            df = df.set_index("id")
    return df


def _featureset_from_features(entity, entity_kwargs):
    """Return featuresets by quering features."""
    results = getattr(query, entity)(**entity_kwargs).all()
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
    suffix: str = None,
    storage_id: str = None,
    created_at=None,
    updated_at=None,
    where: Dict[str, dict] = None,
    as_df: bool = False,
):
    """Query from dobject."""
    model = Table.get_model("dobject")
    kwargs = {k: v for k, v in locals().items() if k in Table.get_fields(model)}
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=model)
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
                dobject_biometas = []
                for biometa in biometas:
                    dobject_biometas += getattr(query, "dobject_biometa")(
                        biometa_id=biometa.id
                    ).all()
                for dobject_biometa in dobject_biometas:
                    dobjects += getattr(query, "dobject")(
                        id=dobject_biometa.dobject_id
                    ).all()
            except AttributeError:
                # query obs metadata
                # find all the link tables to dobject
                dobjects += LinkedQuery().query(
                    entity_return="dobject",
                    entity=entity,
                    entity_kwargs=entity_kwargs,
                )

        results = [i for i in results if i.id in [j.id for j in dobjects]]

    if len(results) > 0:
        for result in results:
            track_usage(result.id, result.v, "query")

    return FilterQueryResultList(model=model, results=results, as_df=as_df)


def query_biometa_from_dobject(
    id: int = None,
    experiment_id: int = None,
    biosample_id: int = None,
    readout_id: int = None,
    featureset_id: int = None,
    dobject_id: str = None,
    as_df: bool = False,
):
    """Query from biometa.

    If dobject_id is provided, will search in the dobject_biometa first.
    """
    model = Table.get_model("biometa")
    kwargs = {k: v for k, v in locals().items() if k in Table.get_fields(model)}
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=model)
    results = _query_stmt(statement=stmt, results_type="all")
    # dobject_id is given, will only return results associated with dobject_id
    if dobject_id is not None:
        biometas = getattr(query, "dobject_biometa")(dobject_id=dobject_id).all()
        if len(biometas) == 0:
            results = biometas
        else:
            biometa_ids = [i.biometa_id for i in biometas]
            results = [i for i in results if i.id in biometa_ids]

    return FilterQueryResultList(model=model, results=results, as_df=as_df)


def _create_query_func(model):
    """Autogenerate query functions for each entity table."""

    def query_func(as_df=False, **kwargs):
        """Query metadata from tables."""
        return Query(model=model, as_df=as_df, kwargs=kwargs)

    query_func.__name__ = model.__name__
    return query_func


class Query:
    def __init__(self, model, as_df=False, kwargs=None) -> None:
        self._model = model
        self._as_df = as_df
        self._kwargs = kwargs

    def _query(self, results_type):
        stmt = _chain_select_stmt(kwargs=self._kwargs, schema_module=self._model)
        results = _query_stmt(statement=stmt, results_type=results_type)
        # track usage for dobjects
        if self._model.__name__ == "dobject":
            for result in results:
                track_usage(result.id, result.v, "query")
        # return DataFrame
        if self._as_df:
            return _return_query_results_as_df(results=results, model=self._model)
        return results

    def all(self):
        return self._query(results_type="all")

    def one(self):
        return self._query(results_type="one")

    def first(self):
        return self._query(results_type="first")


class LinkedQuery:
    """Linked queries."""

    parent_dict = {
        "species": "biosample",
        "biosample_techsample": "biosample",
        "biosample": "biometa",
    }

    def __init__(self) -> None:
        self._engine = settings.instance.db_engine()
        self._inspector = inspect(self._engine)

    @cached_property
    def foreign_keys(self):
        """Foreign keys.

        e.g. {'biosample': {'tissue_id': ('tissue', 'id')}}
        """

        def _get_foreign_keys(table_name, inspector):
            return {
                column["constrained_columns"][0]: (
                    column["referred_table"],
                    column["referred_columns"][0],
                )
                for column in inspector.get_foreign_keys(table_name)
            }

        foreign_keys = {}
        for table_name in self._inspector.get_table_names():
            foreign_keys_table = _get_foreign_keys(table_name, self._inspector)
            if len(foreign_keys_table) > 0:
                foreign_keys[table_name] = foreign_keys_table

        return foreign_keys

    @cached_property
    def foreign_keys_backpop(self):
        """Backpopulated foreign keys.

        e.g. {'tissue': {'id': {'biosample': 'tissue_id'}}}
        """
        foreign_keys_backpop = {}

        for module_name, keys in self.foreign_keys.items():
            for key, (module, ref_key) in keys.items():
                if foreign_keys_backpop.get(module) is None:
                    foreign_keys_backpop[module] = {}
                if foreign_keys_backpop[module].get(ref_key) is None:
                    foreign_keys_backpop[module][ref_key] = {}
                foreign_keys_backpop[module][ref_key][module_name] = key

        return foreign_keys_backpop

    @cached_property
    def link_tables(self):
        """Link tables."""
        link_tables = []
        for name in self._inspector.get_table_names():
            pks = self._inspector.get_pk_constraint(name)["constrained_columns"]
            columns = [i["name"] for i in self._inspector.get_columns(name)]
            if pks == columns and len(self._inspector.get_foreign_keys(name)) > 0:
                link_tables.append(name)

        return link_tables

    def query(self, entity_return, entity, entity_kwargs):
        """Query linked tables via foreign key constraint.

        1. Query fields in entity_n table, whose primary_key is a foreign_key in entity_n-1 table.  # noqa
        2. Query foreign_key in entity_n-1 table, whose primary_key is a foreign_key in entity_n-2 table  # noqa
        3. Repeat until it reaches a linked_table (only contains primary keys).
        """
        results = getattr(query, entity)(**entity_kwargs).all()
        start = entity
        end = entity_return

        current_name = start
        while current_name != end:
            if (
                "id"
                in self._inspector.get_pk_constraint(current_name)[
                    "constrained_columns"
                ]
            ):
                # id is the primary key of current table, aka not a link table
                referred_column = f"{current_name}_id"
                constrained_column = "id"
                parent_name = None
                # if current module id is not present in any other modules as foreign keys  # noqa
                # checks if the any parent module is linked via primary key
                if self.parent_dict.get(current_name) is not None:
                    # specify certain path
                    parent_name = self.parent_dict.get(current_name)
                else:
                    if self.foreign_keys_backpop.get(current_name) is None:
                        for foreign_key in self._inspector.get_foreign_keys(
                            current_name
                        ):
                            if foreign_key["constrained_columns"] == foreign_key[
                                "referred_columns"
                            ] and ("id" in foreign_key["referred_columns"]):
                                parent_name = foreign_key["referred_table"]
                                constrained_column = "id"
                                referred_column = "id"
                    else:
                        parents = self.foreign_keys_backpop.get(current_name)
                        if parents.get("id") is not None:
                            for name, referred_column in parents["id"].items():
                                if referred_column in ["id", "v"]:
                                    continue
                                if name == end:
                                    parent_name = name
                                    break
                                parent_name = name
                parent_results = []
                for result in results:
                    parent_result = getattr(query, parent_name)(
                        **{referred_column: result.__getattribute__(constrained_column)}
                    ).all()
                    parent_results += parent_result
                results = parent_results
            else:
                # if it is a link table to the end module
                if current_name.startswith(f"{end}_"):
                    parent_name = end
                    end_id = f"{end}_id"
                    parent_results = []
                    for result in results:
                        parent_result = getattr(query, end)(
                            **{"id": result.__getattribute__(end_id)}
                        ).all()
                        parent_results += parent_result
                    results = parent_results
                elif self.parent_dict.get(current_name) is not None:
                    parent_name = self.parent_dict.get(current_name)
                    constrained_column = f"{parent_name}_id"
                    parent_results = []
                    for result in results:
                        parent_result = getattr(query, parent_name)(
                            **{"id": result.__getattribute__(constrained_column)}
                        ).all()
                        parent_results += parent_result
                    results = parent_results
                else:
                    pass
            current_name = parent_name

        return results


class FilterQueryResultList:
    def __init__(self, model, results: list, as_df: bool = False) -> None:
        self._model = model
        self._results = results
        self._as_df = as_df

    def _filter(self):
        # track usage for dobjects
        if self._model.__name__ == "dobject":
            for result in self._results:
                track_usage(result.id, result.v, "query")
        # return DataFrame
        if self._as_df:
            return _return_query_results_as_df(results=self._results, model=self._model)
        return self._results

    def all(self):
        return self._filter()

    def one(self):
        if len(self._results) == 0:
            raise exception.NoResultFound
        elif len(self._results) > 1:
            raise exception.MultipleResultsFound
        else:
            results = self._filter()
            if isinstance(results, pd.DataFrame):
                return results.head(1)
            return results[0]

    def first(self):
        if len(self._results) == 0:
            raise exception.NoResultFound
        else:
            results = self._filter()
            if isinstance(results, pd.DataFrame):
                return results.head(1)
            return results[0]


class query:
    """Query literal (semantic) data."""

    pass


for model in Table.list_models():
    func = _create_query_func(model=model)
    setattr(query, model.__name__, staticmethod(func))

setattr(query, "dobject", staticmethod(query_dobject))
setattr(query, "biometa_from_dobject", staticmethod(query_biometa_from_dobject))
