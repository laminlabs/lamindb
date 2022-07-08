from typing import Union
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


def sign_up_hub(user_email) -> Union[str, None]:
    hub = connect_hub()
    secret = id.id_secret()  # generate new secret
    user = hub.auth.sign_up(email=user_email, password=secret)
    # if user already exists a fake user object without identity is returned
    if user.identities:
        logger.info(
            "Please confirm the sign-up email and call again: lndb setup\n\n"
            f"Generated login secret: {secret}.\n"
            "It persists until change of environment or re-install: {settings_file}."
            "Going forward, it is auto-loaded at setup!"
            "You can always recover your secret via your email."
        )
        return secret
    else:
        return None


def sign_in_hub(user_email, secret):
    hub = connect_hub()
    session = hub.auth.sign_in(email=user_email, password=secret)
    data = hub.table("usermeta").select("*").eq("id", session.user.id.hex).execute()
    if len(data.data) > 0:  # user is completely registered
        user_id = data.data[0]["lnid"]
    else:  # user registration on hub gets completed below
        user_id = id.id_user()
        hub.postgrest.auth(session.access_token)
        data = (
            hub.table("usermeta")
            .insert({"id": session.user.id.hex, "lnid": user_id, "handle": user_id})
            .execute()
        )
        assert len(data.data) > 0
    hub.auth.sign_out()
    return user_id
