from datetime import datetime
from enum import Enum
from typing import Optional

from sqlmodel import Field, SQLModel

from ..dev.id import id_file, id_track, id_user


# a user operating the database, e.g., ingesting data
class user(SQLModel, table=True):  # type: ignore
    """Users operating lamindb."""

    id: Optional[str] = Field(default=id_user, primary_key=True)
    name: str
    time_init: datetime = Field(default_factory=datetime.utcnow, nullable=False)


# the process that ingests the data file, the source of the data file
# where the data file comes from
# can be a notebook or a script/pipeline
class interface(SQLModel, table=True):  # type: ignore
    """Interfaces from which `lamindb` is operated."""

    id: str = Field(default=None, primary_key=True)
    name: Optional[str]
    dependency: Optional[str]
    type: str
    user: str = Field(foreign_key="user.id")
    time_init: datetime = Field(default_factory=datetime.utcnow, nullable=False)


# the data file
class file(SQLModel, table=True):  # type: ignore
    """Ingested files storing dense data."""

    id: Optional[str] = Field(default=id_file, primary_key=True)
    name: str
    interface: str = Field(foreign_key="interface.id")
    time_init: datetime = Field(default_factory=datetime.utcnow, nullable=False)


# ----------
# Access log
# ----------


class track_do_type(str, Enum):
    """Data access types."""

    ingest = "ingest"
    query = "query"
    update = "update"
    delete = "delete"
    load = "load"


class track_do(SQLModel, table=True):  # type: ignore
    """Data access log: do operations on the database."""

    id: Optional[str] = Field(default=id_track, primary_key=True)
    type: track_do_type = Field(nullable=False)
    user: str = Field(foreign_key="user.id", nullable=False)
    interface: str = Field(foreign_key="interface.id")
    time: datetime = Field(default_factory=datetime.utcnow, nullable=False)
    file: str = Field(foreign_key="file.id")
