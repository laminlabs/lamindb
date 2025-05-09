from dataclasses import dataclass
from typing import Literal, Optional


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


@dataclass
class TableUID:
    source_table_name: str
    uid_columns: list[str]
    key_constraint: Optional[KeyConstraint]


UIDColumns = list[TableUID]
