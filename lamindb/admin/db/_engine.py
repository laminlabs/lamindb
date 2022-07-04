import sqlmodel as sqm


def get_engine():
    from lamindb.setup import load_settings

    return sqm.create_engine(load_settings().db, future=True)
