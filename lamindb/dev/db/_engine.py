import sqlmodel as sqm


def get_engine(future=True):
    from lamindb._setup import load_or_create_instance_settings

    return sqm.create_engine(load_or_create_instance_settings().db, future=future)
