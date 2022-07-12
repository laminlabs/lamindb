import pandas as pd
from loguru import logger

from lamindb._setup import load_settings
from lamindb._setup._settings import storage_filepath

from .._setup._hub import connect_hub
from ..admin.db import get_engine


def push(dobject_id):
    settings = load_settings()
    # assert settings.cloud_storage

    dobject_metadata = get_dobject_metadata(dobject_id)

    hub = connect_hub()
    session = hub.auth.sign_in(email=settings.user_email, password=settings.user_secret)
    hub.postgrest.auth(session.access_token)

    data = (
        hub.table("dobject")
        .insert({"user_id": session.user.id.hex, **dobject_metadata})
        .execute()
    )
    assert len(data.data) > 0

    logger.info(
        f"Shared dobject ({dobject_id}) from notebook"
        f"by user {settings.user_email} ({settings.user_id})."
    )


def get_dobject_metadata(dobject_id):
    engine = get_engine()

    with engine.connect() as conn:
        df_dobject_metadata = pd.read_sql_table("dobject", conn)
        list_dobject_metadata = df_dobject_metadata[
            df_dobject_metadata["id"] == dobject_id
        ].to_dict("records")
        assert len(list_dobject_metadata) == 1
        dobject_metadata = list_dobject_metadata[0]

    url = storage_filepath(dobject_id)

    return {
        "lnid": dobject_metadata["id"],
        "name": dobject_metadata["name"],
        "suffix": dobject_metadata["suffix"],
        "url": str(url),
    }
