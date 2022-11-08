from typing import Dict, List, Optional, Union

# import bioreadout
# from lamin_logger import colors, logger
from sqlmodel import SQLModel

from lamindb.dev.db._add import add as _add

from .schema._table import table_meta

# from lamindb.dev.db._select import select
# from lamindb.schema import wetlab


def link(
    records1: Union[SQLModel, List[SQLModel]],
    records2: Union[SQLModel, List[SQLModel]],
    *,
    add=True,
) -> Optional[list]:
    """Link data records.

    Create all link records in order to connect entries from two tables.

    Args:
        records1: a single record or a list of records in table1
        records2: a single record or a list of records in table2
        add:
            if True (default): add link records to the db
            if False: return created records without adding to the db

    Returns:
        None or a list of link records
    """
    add_link_records = add  # rename
    if isinstance(records1, SQLModel):
        records1 = [records1]
    if isinstance(records2, SQLModel):
        records2 = [records2]

    # make the two lists the same length
    len1 = len(records1)  # type: ignore
    len2 = len(records2)  # type: ignore
    if len1 == 1:
        records1 = records1 * len2
    if len2 == 1:
        records2 = records2 * len1
    if len(records1) != len(records2):
        raise AssertionError("Can't broadcast the lengths of the two table records!")

    table1_name = records1[0].__table__.name
    table2_name = records2[0].__table__.name

    link_table_name = table_meta.get_link_table(table1_name, table2_name)
    if not link_table_name:
        raise AssertionError(
            f"No link table is found between {table1_name} and {table2_name}!"
        )
    link_table_model = table_meta.get_model(link_table_name)

    # populate the link table
    fks = table_meta.get_foreign_keys(link_table_name)
    fkpks: Dict = {table1_name: [], table2_name: []}
    for k, (f_table, f_id) in fks.items():
        fkpks[f_table].append((k, f_id))

    link_records = []
    for record1, record2 in zip(records1, records2):
        link_record1 = {
            f_id: record1.__getattribute__(k) for f_id, k in fkpks[table1_name]
        }
        link_record2 = {
            f_id: record2.__getattribute__(k) for f_id, k in fkpks[table2_name]
        }

        link_record = link_table_model(**{**link_record1, **link_record2})
        link_records.append(link_record)

    if add_link_records:
        return _add(link_records)
    else:
        return link_records


# def link_readouts(dobject_id: str, efo_id: str):
#     """Link readout to dobjects."""
#     readout = select(wetlab.Readout, efo_id=efo_id).one_or_none()
#     if readout is None:
#         assert sum(i.isdigit() for i in efo_id) == 7
#         readout = add(wetlab.Readout(**bioreadout.readout(efo_id=efo_id)))

#     # select biometa associated with a dobject
#     dobject_biometas = select(wetlab.dobject_biometa, dobject_id=dobject_id).all()
#     if len(dobject_biometas) > 0:
#         biometa_ids = [i.biometa_id for i in dobject_biometas]
#         biometas = (
#             select(wetlab.biometa).where(wetlab.biometa.id.in_(biometa_ids)).all()
#         )
#         for biometa in biometas:
#             biometa.readout_id = readout.id
#         add(biometas)
#     else:
#         add(wetlab.biometa(readout_id=readout.id))
#         logger.warning(
#             f"No biometa found for dobject {dobject_id}, created biometa"
#             f" {biometa_ids[0]}"
#         )

#     logger.success(
#         f"Added {colors.blue(f'readout_id {readout.id}')} to"
#         f" {colors.purple(f'biometa {biometa_ids}')} linked to"
#         f" {colors.green(f'dobject {dobject_id}')}."
#     )
