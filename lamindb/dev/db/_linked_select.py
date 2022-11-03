from functools import cached_property
from typing import Dict

import sqlalchemy as sa
from lndb_setup import settings

from ._core import check_if_link_table, get_foreign_keys
from ._select import select as lnselect


class LinkedSelect:
    """Linked queries."""

    # prefix parent table in cases of confusion
    # the goal here is to reach dobject eventually
    prefix_parents = {
        "bionty.featureset": ["wetlab.biometa"],
        "core.user": ["core.jupynb", "core.pipeline_run"],
        "core.pipeline": ["core.pipeline_run"],
        "core.pipeline_run": ["core.dtransform"],
        "core.dtransform": ["core.dobject"],
        "core.dtransform_in": ["core.dtransform"],
        "wetlab.biosample_techsample": ["wetlab.biosample"],
        "bionty.species": ["wetlab.biosample"],
        "wetlab.biosample": ["wetlab.biometa"],
    }

    def __init__(self) -> None:
        self._inspector = sa.inspect(settings.instance.db_engine())

    @cached_property
    def foreign_keys_by_tables(self):
        """Foreign keys for all tables indexed by tables.

        Foreign keys are parents both in two ways:

        1. One queries a table by foreign key to select a value for a primary key.
        2. One typically needs the foreign key first to insert a row in a table
           that refers to it.

        Returns: {
            'dobject_biometa':
                {'dobject_id': ('dobject', 'id'), 'biometa_id': ('biometa', 'id')}
        }
        """
        foreign_keys = {}
        for table_name in self._inspector.get_table_names():
            foreign_keys_dict = get_foreign_keys(table_name, self._inspector)
            if len(foreign_keys_dict) > 0:
                foreign_keys[table_name] = foreign_keys_dict
        return foreign_keys

    @cached_property
    def tables_by_foreign_keys(self):
        """All tables for all foreign keys indexed by foreign keys.

        Returns: {
            'user':
                {'id': [('jupynb', 'created_by'), ('pipeline_run', 'created_by')]
        }
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
            return None

        pkfks = check_if_link_table(table_name)
        if pkfks:
            # link tables, meaning the parent table shares a primary key
            fks = {
                k: [v]
                for k, v in self.foreign_keys_by_tables.get(table_name).items()
                if k in pkfks
            }
        else:
            # tables that are linked to the parent table via foreign key constraints
            # these are the cases which id is a column called `{table}_id` in the parent table  # noqa
            fks = self.tables_by_foreign_keys.get(table_name)

        if fks:
            results: Dict = {}
            for cons_col, parents in fks.items():
                for name, ref_col in parents:
                    if name not in results:
                        results[name] = []
                    results[name].append((cons_col, ref_col))

            # we prefix certain table's parents when there is > 1 parents
            prefix_parents = self.prefix_parents.get(table_name)
            if table_name.startswith("featureset_"):
                prefix_parents = ["featureset"]
            if table_name.startswith("dobject_"):
                prefix_parents = ["dobject"]
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
                parent_results[parent_name] += lnselect(
                    parent_name,
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
        results = lnselect(table_name, **table_kwargs).all()
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
