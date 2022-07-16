from datetime import datetime as datetime
from enum import Enum
from typing import Optional

from sqlmodel import Field, ForeignKeyConstraint, SQLModel, UniqueConstraint

from ..dev.id import id_dobject, id_track


def utcnow():
    return datetime.utcnow().replace(microsecond=0)


class user(SQLModel, table=True):  # type: ignore
    """Users operating `lamindb`."""

    __table_args__ = (UniqueConstraint("email"),)
    id: Optional[str] = Field(primary_key=True)
    email: str
    time_created: datetime = Field(default_factory=utcnow, nullable=False)
    time_updated: datetime = Field(default_factory=utcnow, nullable=False)


class interface(SQLModel, table=True):  # type: ignore
    """Interface from which users & machines operate `lamindb`."""

    id: str = Field(default=None, primary_key=True)
    v: str = Field(default=None, primary_key=True)
    name: Optional[str]
    type: str  #: Jupyter notebook (nbproject), pipeline, etc.
    user_id: str = Field(foreign_key="user.id")
    time_created: datetime = Field(default_factory=utcnow, nullable=False)
    time_updated: datetime = Field(default_factory=utcnow, nullable=False)


class dobject(SQLModel, table=True):  # type: ignore
    """Data objects in storage & memory.

    Data objects often correspond to files.

    Storage ⟷ memory examples:
    - `.csv`, `.tsv`, `.feather`, `.parquet` ⟷ `pd.DataFrame`
    - `.h5ad`, `.h5mu`, or their zarr versions ⟷ `anndata.AnnData`, `mudata.MuData`
    - `.jpg`, `.png` ⟷ `np.ndarray`, or a dedicated imaging in-memory container
    - zarr directory ⟷ zarr loader
    - TileDB store ⟷ TileDB loader
    - fastq ⟷ ?
    - .vcf ⟷ ?
    """

    __table_args__ = (
        ForeignKeyConstraint(
            ["interface_id", "interface_v"],
            ["interface.id", "interface.v"],
        ),
    )

    id: Optional[str] = Field(default_factory=id_dobject, primary_key=True)
    v: str = Field(default=None, primary_key=True)
    name: Optional[str]
    file_suffix: str
    interface_id: str = Field(foreign_key="interface.id")
    interface_v: str = Field(foreign_key="interface.v")
    time_created: datetime = Field(default_factory=utcnow, nullable=False)
    time_updated: datetime = Field(default_factory=utcnow, nullable=False)


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
