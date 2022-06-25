import sqlmodel as sqm


def get_engine():
    from lamindb.setup import settings

    return sqm.create_engine(settings().db, future=True)
