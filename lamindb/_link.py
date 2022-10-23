from typing import Dict, List, Optional, Union

from sqlmodel import SQLModel

from .dev.db._add import add
from .dev.db._core import get_foreign_keys, get_link_table
from .schema._table import Table


def link(
    entries_table1: Union[SQLModel, List[SQLModel]],
    entries_table2: Union[SQLModel, List[SQLModel]],
    *,
    add_link_entries=True,
) -> Optional[list]:
    """Populate link table entries of two entities.

    Args:
        entries_table1: a single entry or a list of entries in table1
        entries_table2: a single entry or a list of entries in table2
        add_link_entries:
            if True (default): add link entries to the db
            if False: return created entries without adding to the db

    Returns:
        None or a list of link entries
    """
    if isinstance(entries_table1, SQLModel):
        entries_table1 = [entries_table1]
    if isinstance(entries_table2, SQLModel):
        entries_table2 = [entries_table2]

    # make the two lists the same length
    len1 = len(entries_table1)  # type: ignore
    len2 = len(entries_table2)  # type: ignore
    if len1 == 1:
        entries_table1 = entries_table1 * len2
    if len2 == 1:
        entries_table2 = entries_table2 * len1
    if len(entries_table1) != len(entries_table2):
        raise AssertionError("Can't broadcast the lengths of the two table entries!")

    table1_name = entries_table1[0].__table__.name
    table2_name = entries_table2[0].__table__.name

    link_table_name = get_link_table(table1_name, table2_name)
    if not link_table_name:
        raise AssertionError(
            f"No link table is found between {table1_name} and {table2_name}!"
        )
    link_table_model = Table.get_model(link_table_name)

    # populate the link table
    fks = get_foreign_keys(link_table_name)
    fkpks: Dict = {table1_name: [], table2_name: []}
    for k, (f_table, f_id) in fks.items():
        fkpks[f_table].append((k, f_id))

    link_entries = []
    for entry1, entry2 in zip(entries_table1, entries_table2):
        link_entry1 = {
            f_id: entry1.__getattribute__(k) for f_id, k in fkpks[table1_name]
        }
        link_entry2 = {
            f_id: entry2.__getattribute__(k) for f_id, k in fkpks[table2_name]
        }

        link_entry = link_table_model(**{**link_entry1, **link_entry2})
        link_entries.append(link_entry)
        # TODO: remove after the new add API
        kwargs = link_entry.dict()
        if add_link_entries:
            add(link_table_model(**kwargs))

    if add_link_entries:
        return None
    else:
        return link_entries
