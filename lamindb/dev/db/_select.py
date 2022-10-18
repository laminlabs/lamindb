from datetime import datetime
from functools import cached_property
from typing import Dict, Optional

import bionty as bt
import sqlalchemy as sa
import sqlmodel as sqm
from lndb_setup import settings

from ...schema._table import Table
from ._core import check_if_link_table, get_foreign_keys
from ._select_result import SelectResult
from ._track_usage import track_usage


def _select_stmt(statement, results_type="all"):
    with sqm.Session(settings.instance.db_engine()) as session:
        results = session.exec(statement).__getattribute__(results_type)()
    return results


def _chain_select_stmt(kwargs: dict, schema_module):
    stmt = sqm.select(schema_module)
    if len(kwargs) > 0:
        for field, value in kwargs.items():
            if (field in {"cls", "self"}) or (value is None):
                continue
            else:
                stmt = stmt.where(getattr(schema_module, field) == value)
    return stmt


def _featureset_from_features(entity, entity_kwargs):
    """Return featuresets by quering features."""
    results = getattr(select, entity)(**entity_kwargs).all()
    if len(results) > 0:
        featureset_ids = []
        for feature in results:
            featuresets = getattr(select, f"featureset_{entity}")(
                **{f"{entity}_id": feature.id}
            ).all()
            if len(featuresets) > 0:
                featureset_ids += [
                    featureset.featureset_id for featureset in featuresets
                ]
        return list(set(featureset_ids))
    return []


def _create_select_func(model):
    """Autogenerate select functions for each entity table."""

    def select_func(**kwargs):
        return _select(model=model, kwargs=kwargs)

    select_func.__name__ = model.__name__
    import_module, prefix = (
        model.__module__.split(".")[0],
        model.__module__.split(".")[1],
    )
    prefix = f".{prefix.replace('_schema', 'schema')}" if "schema" in prefix else ""
    url = f"https://lamin.ai/docs/{import_module.replace('_', '-')}/{import_module}{prefix}.{model.__name__}"  # noqa
    select_func.__doc__ = (
        f"""Select metadata from `{import_module}.{model.__name__} <{url}>`__."""
    )
    return select_func


def _select(model, result_list=None, kwargs=dict()) -> SelectResult:
    """Simple queries."""
    if result_list is not None:
        return SelectResult(results=result_list, model=model)

    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=model)
    results = _select_stmt(statement=stmt, results_type="all")
    # track usage for dobjects
    if model.__name__ == "dobject":
        for result in results:
            track_usage(result.id, "select")

    return SelectResult(results=results, model=model)


def select_dobject(
    id: Optional[str] = None,
    v: Optional[str] = None,
    name: Optional[str] = None,
    size: Optional[float] = None,
    dtransform_id: Optional[str] = None,
    suffix: Optional[str] = None,
    storage_id: Optional[str] = None,
    checksum: Optional[str] = None,
    created_at: datetime = None,
    updated_at: datetime = None,
    where: Dict[str, dict] = None,
):
    """Select from dobject.

    `lnschema_core.dobject <https://lamin.ai/docs/lnschema-core/lnschema_core.dobject>`__.
    """  # noqa
    model = Table.get_model("dobject")
    kwargs = {k: v for k, v in locals().items() if k in Table.get_fields(model)}
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=model)
    results = _select_stmt(statement=stmt, results_type="all")

    if where is not None:
        dobjects = []
        for entity, entity_kwargs in where.items():
            # TODO: this part needs refactor
            try:
                bt.lookup.feature_model.__getattribute__(entity)
                # select features
                featureset_ids = _featureset_from_features(
                    entity=entity, entity_kwargs=entity_kwargs
                )
                biometas = []
                for featureset_id in featureset_ids:
                    biometas += getattr(select, "biometa")(
                        featureset_id=featureset_id
                    ).all()
                dobject_biometas = []
                for biometa in biometas:
                    dobject_biometas += getattr(select, "dobject_biometa")(
                        biometa_id=biometa.id
                    ).all()
                for dobject_biometa in dobject_biometas:
                    dobjects += getattr(select, "dobject")(
                        id=dobject_biometa.dobject_id
                    ).all()
            except AttributeError:
                # select obs metadata
                # find all the link tables to dobject
                dobjects += LinkedSelect().select(
                    table_name=entity,
                    table_kwargs=entity_kwargs,
                    return_table="dobject",
                )

        results = [i for i in results if i.id in [j.id for j in dobjects]]

    return _select(model=model, result_list=results)


class LinkedSelect:
    """Linked queries."""

    # prefixed parent table in cases of confusion
    # the goal here is to reach dobject eventually
    prefix_parents = {
        "featureset": ["biometa"],
        "user": ["jupynb", "pipeline_run"],
        "pipeline": ["pipeline_run"],
        "featureset_cell_marker": ["featureset"],
        "featureset_gene": ["featureset"],
        "featureset_protein": ["featureset"],
        "pipeline_run": ["dtransform"],
        "dtransform": ["dobject"],
        "dobject_bfxmeta": ["dobject"],
        "dobject_biometa": ["dobject"],
        "dtransform_in": ["dtransform"],
        "biosample_techsample": ["biosample"],
        "species": ["biosample"],
        "biosample": ["biometa"],
    }

    def __init__(self) -> None:
        self._engine = settings.instance.db_engine()
        self._inspector = sa.inspect(self._engine)

    @cached_property
    def foreign_keys(self):
        """Foreign keys of all tables.

        Returns: {'dobject_biometa':
                {'dobject_id': ('dobject', 'id'), 'biometa_id': ('biometa', 'id')}}
        """
        foreign_keys = {}
        for table_name in self._inspector.get_table_names():
            foreign_keys_dict = get_foreign_keys(table_name, self._inspector)
            if len(foreign_keys_dict) > 0:
                foreign_keys[table_name] = foreign_keys_dict

        return foreign_keys

    @cached_property
    def foreign_keys_backpop(self):
        """Backpopulated foreign keys of all tables.

        Returns: {'user':
                {'id': [('jupynb', 'created_by'), ('pipeline_run', 'created_by')]}
        """
        results = {}

        for table, keys in self.foreign_keys.items():
            for cons_key, (ref_table, ref_key) in keys.items():
                if results.get(ref_table) is None:
                    results[ref_table] = {}
                if results[ref_table].get(ref_key) is None:
                    results[ref_table][ref_key] = []
                results[ref_table][ref_key].append((table, cons_key))

        return results

    def get_parent_tables(self, table_name: str):
        """Return all tables containing a column with foreign key to the provided table.

        Returns {parent_table_name : [(constraint_column, referred_column)]}
        """
        # skip migration and version tables
        # dobject should be the end table
        if table_name.startswith(("migration_", "version_")) or table_name in [
            "dobject"
        ]:
            return

        pkfks = check_if_link_table(table_name)
        if pkfks:
            # link tables, meaning the parent table shares a primary key
            fks = {
                k: [v]
                for k, v in self.foreign_keys.get(table_name).items()
                if k in pkfks
            }
        else:
            # tables that are linked to the parent table via foreign key constraints
            # these are the cases which id is a column called `{table}_id` in the parent table  # noqa
            fks = self.foreign_keys_backpop.get(table_name)

        if fks:
            results: Dict = {}
            for cons_col, parents in fks.items():
                for name, ref_col in parents:
                    if name not in results:
                        results[name] = []
                    results[name].append((cons_col, ref_col))

            # we prefix certain table's parents when there is > 1 parents
            prefix_parents = self.prefix_parents.get(table_name)
            if len(results) > 1 and prefix_parents:
                return {k: v for k, v in results.items() if k in prefix_parents}

            return results

    def select_from_parents(self, results: list, table_name: str) -> list:
        """Returns select results from parent tables."""
        if not results:
            raise AssertionError("No results found!")
        parents = self.get_parent_tables(table_name)

        parent_results: Dict = {}
        for parent_name, keys in parents.items():
            if parent_name not in parent_results:
                parent_results[parent_name] = []
            for result in results:
                parent_results[parent_name] += getattr(select, parent_name)(
                    **{
                        ref_col: result.__getattribute__(const_col)
                        for const_col, ref_col in keys
                    }
                ).all()

        return [(k, results) for k, results in parent_results.items()]

    def select_from_single_parent(self, results: list, table_name: str, end_table: str):
        current_table = table_name
        while current_table != end_table:
            parent_results_all = self.select_from_parents(
                results=results, table_name=current_table
            )
            # single parent
            if len(parent_results_all) == 1:
                parent_name, parent_results = parent_results_all[0]
                current_table = parent_name
                results = parent_results
            if len(parent_results_all) > 1:
                return {k: results for k, results in parent_results_all}
        return results

    def get_pks(self, table_name: str) -> list:
        """Return a list of primary keys."""
        return self._inspector.get_pk_constraint(table_name)

    def select(self, table_name: str, table_kwargs: dict, return_table: str):
        """Select from linked tables."""
        results = getattr(select, table_name)(**table_kwargs).all()
        linked_results = self.select_from_single_parent(
            results=results, table_name=table_name, end_table=return_table
        )
        if isinstance(linked_results, dict):
            combine_linked_results = []
            for table_name, results in linked_results.items():
                if not results:
                    continue
                results_tmp = self.select_from_single_parent(
                    results=results,
                    table_name=table_name,
                    end_table=return_table,
                )
                combine_linked_results += [
                    i for i in results_tmp if i not in combine_linked_results
                ]

            return combine_linked_results
        return list(linked_results)


class select:
    """Select data.

    Data is queried by entity.
    Arguments to entity allow to express deep linked queries very concisely.

    Guide: :doc:`/db/guide/select-load`.

    Returns a :class:`~lamindb.dev.db.SelectResult` object.
    """

    pass


for model in Table.list_models():
    func = _create_select_func(model=model)
    setattr(select, model.__name__, staticmethod(func))

setattr(select, "dobject", staticmethod(select_dobject))
