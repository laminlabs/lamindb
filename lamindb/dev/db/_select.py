from datetime import datetime
from functools import cached_property
from typing import Dict, Optional

import bionty as bt
import sqlalchemy as sa
import sqlmodel as sqm
from lndb_setup import settings

from ...schema._table import Table
from ._core import get_foreign_keys, get_link_tables
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
                    entity_return="dobject",
                    entity=entity,
                    entity_kwargs=entity_kwargs,
                )

        results = [i for i in results if i.id in [j.id for j in dobjects]]

    return _select(model=model, result_list=results)


class LinkedSelect:
    """Linked queries."""

    parent_dict = {
        "species": "biosample",
        "biosample_techsample": "biosample",
        "biosample": "biometa",
    }

    def __init__(self) -> None:
        self._engine = settings.instance.db_engine()
        self._inspector = sa.inspect(self._engine)

    @cached_property
    def foreign_keys(self):
        """Foreign keys.

        e.g. {'biosample': {'tissue_id': ('tissue', 'id')}}
        """
        foreign_keys = {}
        for table_name in self._inspector.get_table_names():
            foreign_keys_table = get_foreign_keys(table_name, self._inspector)
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
        return get_link_tables(self._inspector)

    def select(self, entity_return, entity, entity_kwargs):
        """Select linked tables via foreign key constraint.

        1. Select fields in entity_n table, whose primary_key is a foreign_key in entity_n-1 table.  # noqa
        2. Select foreign_key in entity_n-1 table, whose primary_key is a foreign_key in entity_n-2 table  # noqa
        3. Repeat until it reaches a linked_table (only contains primary keys).
        """
        results = getattr(select, entity)(**entity_kwargs).all()
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
                    parent_result = getattr(select, parent_name)(
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
                        parent_result = getattr(select, end)(
                            **{"id": result.__getattribute__(end_id)}
                        ).all()
                        parent_results += parent_result
                    results = parent_results
                elif self.parent_dict.get(current_name) is not None:
                    parent_name = self.parent_dict.get(current_name)
                    constrained_column = f"{parent_name}_id"
                    parent_results = []
                    for result in results:
                        parent_result = getattr(select, parent_name)(
                            **{"id": result.__getattribute__(constrained_column)}
                        ).all()
                        parent_results += parent_result
                    results = parent_results
                else:
                    pass
            current_name = parent_name

        return results


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
