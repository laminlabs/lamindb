from typing import Dict, List, Optional, Union

import bioreadout
from lamin_logger import colors, logger
from sqlmodel import SQLModel

from lamindb.dev.db._add import add
from lamindb.dev.db._select import select
from lamindb.schema import wetlab

from .dev.db._core import get_foreign_keys, get_link_table
from .schema._table import table_meta


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
    link_table_model = table_meta.get_model(link_table_name)

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


def link_readouts(dobject_id: str, efo_id: str):
    """Link readout to dobjects."""
    readout = select(wetlab.readout, efo_id=efo_id).one_or_none()
    if readout is None:
        assert sum(i.isdigit() for i in efo_id) == 7
        readout = add(wetlab.readout(**bioreadout.readout(efo_id=efo_id)))

    # select biometa associated with a dobject
    dobject_biometas = select(wetlab.dobject_biometa, dobject_id=dobject_id).all()
    if len(dobject_biometas) > 0:
        biometa_ids = [i.biometa_id for i in dobject_biometas]
        biometas = (
            select(wetlab.biometa).where(wetlab.biometa.id.in_(biometa_ids)).all()
        )
        for biometa in biometas:
            biometa.readout_id = readout.id
        add(biometas)
    else:
        add(wetlab.biometa(readout_id=readout.id))
        logger.warning(
            f"No biometa found for dobject {dobject_id}, created biometa"
            f" {biometa_ids[0]}"
        )

    logger.success(
        f"Added {colors.blue(f'readout_id {readout.id}')} to"
        f" {colors.purple(f'biometa {biometa_ids}')} linked to"
        f" {colors.green(f'dobject {dobject_id}')}."
    )
