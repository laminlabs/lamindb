from typing import List, Optional
from uuid import uuid4

from sqlmodel import Field, Relationship, SQLModel


def uuid4_char32():
    """Convert to hex string as Python UUID type has no stable conversion to SQL."""
    return uuid4().hex


class Publication(SQLModel, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True, nullable=False)
    uuid: Optional[str] = Field(
        default_factory=uuid4_char32,
        index=True,
        nullable=False,
    )
    name: str = Field(index=True)
    dataset: List["Dataset"] = Relationship(back_populates="publication")


class Dataset(SQLModel, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True, nullable=False)
    uuid: Optional[str] = Field(
        default_factory=uuid4_char32,
        index=True,
        nullable=False,
    )
    name: str = Field(index=True)
    publication_id: Optional[int] = Field(default=None, foreign_key="publication.id")
    publication: Optional[Publication] = Relationship(back_populates="dataset")
