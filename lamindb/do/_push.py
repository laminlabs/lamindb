import uuid

from lamin_logger import logger
from lndb_setup import settings
from lndb_setup._hub import connect_hub
from supabase import Client

from ._load import load


def get_hub_with_authentication():
    hub = connect_hub()
    session = hub.auth.sign_in(
        email=settings.user.user_email, password=settings.user.user_secret
    )
    hub.postgrest.auth(session.access_token)
    return hub


class hub:
    """Access the hub."""

    @classmethod
    def share_instance(cls):
        """Publish instance with all dobjects."""
        hub = get_hub_with_authentication()
        dobject_df = load("dobject")
        for id, v in dobject_df.index:
            cls.share_dobject(id, v, hub)

    @classmethod
    def share_dobject(cls, id: str, v: str, hub: Client = None):
        """Publish a single dobject."""
        if hub is None:
            hub = get_hub_with_authentication()
        instance = get_or_create_instance(hub)
        dobject = load("dobject").loc[[(id, v)]].reset_index().to_dict("records")[0]
        if not jupynb_exists(dobject["jupynb_id"], dobject["jupynb_v"], hub):
            logger.info(
                f'Published jupynb ({dobject["jupynb_id"]}, {dobject["jupynb_v"]}).'
            )
            insert_jupynb(
                dobject["jupynb_id"], dobject["jupynb_v"], instance["id"], hub
            )
        if not dobject_exists(id, v, hub):
            logger.info(f"Published dobject ({id}, {v}).")
            insert_dobject(id, v, hub)
        if not user_dobject_exists(id, v, hub):
            insert_user_dobject(id, v, hub)

    @classmethod
    def delete_instance(cls):
        """Delete the instance on the hub."""
        hub = get_hub_with_authentication()
        instance = get_or_create_instance(hub)
        dobject_df = load("dobject")
        jupynb_df = load("jupynb")
        for id, v in dobject_df.index:
            logger.info(f"Deleted dobject ({id}, {v}).")
            delete_user_dobject(id, v, hub)
            delete_dobject(id, v, hub)
        for id, v in jupynb_df.index:
            logger.info(f"Deleted jupynb ({id}, {v}).")
            delete_jupynb(id, v, hub)
        delete_user_instance(instance["id"], hub)
        delete_instance(hub)


# Instance


def get_or_create_instance(hub: Client):
    instance = get_instance(hub)
    if instance is None:
        instance = insert_instance(hub)
        logger.info("Publishing instance.")
    if not user_instance_exists(instance["id"], hub):
        insert_user_instance(instance["id"], hub)
    return instance


def get_instance(hub: Client):
    data = (
        hub.table("instance")
        .select("*")
        .eq("name", settings.instance.instance_name)
        .execute()
    )
    if len(data.data) > 0:
        return data.data[0]
    return None


def delete_instance(hub: Client):
    (
        hub.table("instance")
        .delete()
        .eq("name", settings.instance.instance_name)
        .execute()
    )


def insert_instance(hub: Client):
    data = (
        hub.table("instance")
        .insert(
            {
                "id": str(uuid.uuid4()),
                "name": settings.instance.instance_name,
                "owner_id": hub.auth.session().user.id.hex,
                "storage_dir": str(settings.instance.storage_dir),
                "dbconfig": settings.instance._dbconfig,
                "cache_dir": str(settings.instance.cache_dir),
                "sqlite_file": str(settings.instance._sqlite_file),
                "sqlite_file_local": str(settings.instance._sqlite_file_local),
                "db": settings.instance.db,
            }
        )
        .execute()
    )
    assert len(data.data) == 1
    return data.data[0]


# User instance


def user_instance_exists(instance_id, hub: Client):
    data = (
        hub.table("user_instance")
        .select("*")
        .eq("user_id", hub.auth.session().user.id.hex)
        .eq("instance_id", instance_id)
        .execute()
    )
    if len(data.data) > 0:
        return True
    return False


def delete_user_instance(instance_id, hub: Client):
    (
        hub.table("user_instance")
        .delete()
        .eq("user_id", hub.auth.session().user.id.hex)
        .eq("instance_id", instance_id)
        .execute()
    )


def insert_user_instance(instance_id, hub: Client):
    data = (
        hub.table("user_instance")
        .insert(
            {
                "user_id": hub.auth.session().user.id.hex,
                "instance_id": instance_id,
            }
        )
        .execute()
    )
    assert len(data.data) == 1
    return data.data[0]


# jupynb


def jupynb_exists(id, v, hub: Client):
    data = hub.table("jupynb").select("*").eq("id", id).eq("v", v).execute()
    if len(data.data) > 0:
        return True
    return False


def delete_jupynb(id, v, hub: Client):
    hub.table("jupynb").delete().eq("id", id).eq("v", v).execute()


def insert_jupynb(id, v, instance_id, hub: Client):
    jupynb = load("jupynb").loc[[(id, v)]].reset_index().to_dict("records")[0]
    data = (
        hub.table("jupynb")
        .insert(
            {
                "id": jupynb["id"],
                "v": jupynb["v"],
                "name": jupynb["name"],
                "type": jupynb["type"],
                "instance_id": instance_id,
                "owner_id": hub.auth.session().user.id.hex,
                # We could use the get_user function the hub user id from the
                # lnid of the original creator of a notebook, but this won't
                # works because RLS policy. So for the moment the owner_id is
                # the id of the first user that push a dobject from this
                # notebook.
                "time_created": jupynb["time_created"].strftime("%Y-%m-%d %H:%M:%S"),
                "time_updated": jupynb["time_updated"].strftime("%Y-%m-%d %H:%M:%S"),
            }
        )
        .execute()
    )
    assert len(data.data) == 1
    return data.data[0]


# dobject


def dobject_exists(id, v, hub: Client):
    data = hub.table("dobject").select("*").eq("id", id).eq("v", v).execute()
    if len(data.data) > 0:
        return True
    return False


def delete_dobject(id, v, hub: Client):
    hub.table("dobject").delete().eq("id", id).eq("v", v).execute()


def insert_dobject(id, v, hub: Client):
    dobject = load("dobject").loc[[(id, v)]].reset_index().to_dict("records")[0]
    data = (
        hub.table("dobject")
        .insert(
            {
                "id": dobject["id"],
                "v": dobject["v"],
                "name": dobject["name"],
                "file_suffix": dobject["file_suffix"],
                "jupynb_id": dobject["jupynb_id"],
                "jupynb_v": dobject["jupynb_v"],
                "time_created": dobject["time_created"].strftime("%Y-%m-%d %H:%M:%S"),
                "time_updated": dobject["time_updated"].strftime("%Y-%m-%d %H:%M:%S"),
            }
        )
        .execute()
    )
    assert len(data.data) == 1
    return data.data[0]


# User dobject


def user_dobject_exists(id, v, hub: Client):
    data = (
        hub.table("user_dobject")
        .select("*")
        .eq("dobject_id", id)
        .eq("dobject_v", v)
        .eq("user_id", hub.auth.session().user.id.hex)
        .execute()
    )
    if len(data.data) > 0:
        return True
    return False


def delete_user_dobject(id, v, hub: Client):
    (
        hub.table("user_dobject")
        .delete()
        .eq("dobject_id", id)
        .eq("dobject_v", v)
        .eq("user_id", hub.auth.session().user.id.hex)
        .execute()
    )


def insert_user_dobject(id, v, hub: Client):
    data = (
        hub.table("user_dobject")
        .insert(
            {
                "dobject_id": id,
                "dobject_v": v,
                "user_id": hub.auth.session().user.id.hex,
            }
        )
        .execute()
    )
    assert len(data.data) == 1
    return data.data[0]
