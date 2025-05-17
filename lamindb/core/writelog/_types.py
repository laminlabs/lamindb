from dataclasses import dataclass
from enum import Enum
from typing import Literal, Optional


class ColumnType(Enum):
    INT = 0
    BOOL = 1
    STR = 2
    DATE = 3
    FLOAT = 4
    JSON = 5
    TIMESTAMPTZ = 6


@dataclass(frozen=True)
class Column:
    name: str
    type: ColumnType
    ordinal_position: int


@dataclass
class KeyConstraint:
    """Simple encapsulation of one of a table's key constraints."""

    constraint_name: str
    constraint_type: Literal["PRIMARY KEY", "FOREIGN KEY"]

    # These need to be a list to account for composite primary keys
    source_columns: list[Column]
    target_columns: list[Column]

    target_table: str


@dataclass
class TableUID:
    source_table_name: str
    uid_columns: list[str]
    key_constraint: Optional[KeyConstraint]


UIDColumns = list[TableUID]
