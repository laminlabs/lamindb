import pandas as pd
from loguru import logger
from supabase import create_client

from lamindb._setup import load_settings
from lamindb._setup._settings import storage_filepath

from ..admin.db import get_engine


def push(dobject_id):
    settings = load_settings()
    # assert settings.cloud_storage

    dobject_metadata = get_dobject_metadata(dobject_id)

    hub = create_client(
        "https://qntuvxhregqtypdorhaw.supabase.co",
        "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9"
        ".eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InFudHV2eGhy"
        "ZWdxdHlwZG9yaGF3Iiwicm9sZSI6ImFub24iLCJpYXQiOj"
        "E2NTcyMjMxODksImV4cCI6MTk3Mjc5OTE4OX0"
        ".hzrRK-xbhFSqFb1-cVu6NFM-UcPON4HBdaH2qe3vKbA",
    )
    session = hub.auth.sign_in(
        email="frederic.enard@gmail.com", password=settings.user_secret
    )
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
        "name": dobject_metadata["name"],
        "suffix": dobject_metadata["suffix"],
        "url": str(url),
    }
