from urllib.request import urlretrieve

from supabase import create_client

from .._logger import logger
from ..dev import id
from ._settings_store import Connector


def connect_hub():
    connector_file, _ = urlretrieve(
        "https://lamin-site-assets.s3.amazonaws.com/connector.env"
    )
    connector = Connector(_env_file=connector_file)
    return create_client(connector.url, connector.key)


def signup_hub(user_email):
    hub = connect_hub()
    secret = id.id_secret()
    hub.auth.sign_up(email=user_email, password=secret)
    logger.info(
        "Please confirm the sign-up email and then repeat the login call.\n"
        f"Generated login secret: {secret}.\n"
        "Your secret has been stored and will persist until a new installation."
    )
    return secret


def signin_hub(user_email, secret=None):
    hub = connect_hub()
    session = hub.auth.sign_in(email=user_email, password=secret)
    data = hub.table("usermeta").select("*").eq("id", session.user.id.hex).execute()
    if len(data.data) > 0:
        user_id = data.data[0]["lnid"]
    else:
        user_id = id.id_user()
        hub.postgrest.auth(session.access_token)
        data = (
            hub.table("usermeta")
            .insert({"id": session.user.id.hex, "lnid": user_id, "handle": user_id})
            .execute()
        )
        assert len(data.data) > 0

    hub.auth.sign_out()
