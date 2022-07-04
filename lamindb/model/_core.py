from datetime import datetime as datetime
from enum import Enum
from typing import Optional

from sqlmodel import Field, SQLModel

from ..dev.id import id_dobject, id_track, id_user


def utcnow():
    return datetime.utcnow().replace(microsecond=0)


class user(SQLModel, table=True):  # type: ignore
    """Users operating `lamindb`."""

    id: Optional[str] = Field(default_factory=id_user, primary_key=True)
    name: str
    time_init: datetime = Field(default_factory=utcnow, nullable=False)


class interface(SQLModel, table=True):  # type: ignore
    """User interfaces from which users operate `lamindb`."""

    id: str = Field(default=None, primary_key=True)
    name: Optional[str]
    dependency: Optional[str]  #: Environment dependencies of operation.
    type: str  #: Jupyter notebook, pipelines, etc.
    user_id: str = Field(foreign_key="user.id")
    time_init: datetime = Field(default_factory=utcnow, nullable=False)


class dobject(SQLModel, table=True):  # type: ignore
    """Data objects in storage & memory, often correspond to files.

    Storage ⟷ memory examples:
    - `.csv`, `.tsv`, `.feather`, `.parquet` ⟷ `pd.DataFrame`
    - `.h5ad`, `.h5mu` ⟷ `anndata.AnnData`, `mudata.MuData`
    - `.jpg`, `.png` ⟷ `np.ndarray`
    - Zarr directory ⟷ Zarr loader
    - TileDB store ⟷ TileDB loader
    """

    id: Optional[str] = Field(default_factory=id_dobject, primary_key=True)
    name: str
    interface_id: str = Field(foreign_key="interface.id")
    time_init: datetime = Field(default_factory=utcnow, nullable=False)


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

    id: Optional[str] = Field(default_factory=id_track, primary_key=True)
    type: track_do_type = Field(nullable=False)  #: Data access type.
    user_id: str = Field(foreign_key="user.id", nullable=False)
    interface_id: str = Field(foreign_key="interface.id")
    time: datetime = Field(default_factory=utcnow, nullable=False)
    dobject_id: str = Field(foreign_key="dobject.id")
