import sqlmodel as sqm

from lamindb.setup import settings


def get_engine():
    return sqm.create_engine(settings().db, future=True)
