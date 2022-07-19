from urllib.request import urlretrieve

from pydantic import BaseSettings


def get_default_supabase_credentials():
    connector_file, _ = urlretrieve(
        "https://lamin-site-assets.s3.amazonaws.com/connector.env"
    )
    connector = Connector(_env_file=connector_file)
    return connector.url, connector.key


class Connector(BaseSettings):
    url: str
    key: str
