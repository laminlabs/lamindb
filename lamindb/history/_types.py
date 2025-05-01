from dataclasses import dataclass
from typing import Literal


@dataclass
class ManyToManyRelationship:
    first_column: str
    first_table: str
    second_column: str
    second_table: str


@dataclass
class KeyConstraint:
    """Simple encapsulation of one of a table's key constraints."""

    constraint_name: str
    constraint_type: Literal["PRIMARY KEY", "FOREIGN KEY"]

    # These need to be a list to account for composite primary keys
    source_columns: list[str]
    target_columns: list[str]

    target_table: str


# Record UIDs are a mapping from table name to a list of column names.
# When we actually store history, we'll store table IDs instead of names.
# For most tables, record UIDs will look like {<table's ID>: [list of table's
# UID columns]}. For many-to-many tables, the record's UID will be a mapping
# of the linked tables' IDs to their respective UID columns.
UIDColumns = dict[str, list[str]]
