import sqlmodel as sqm


def get_engine():
    from lamindb._setup import load_instance_settings

    return sqm.create_engine(load_instance_settings().db, future=True)
