import uuid

from lndb_setup._hub import connect_hub
from lndb_setup._setup_instance import (
    load_or_create_instance_settings,
    load_or_create_user_settings,
)

from ._load import load


def get_hub_with_authentication():
    user_settings = load_or_create_user_settings()
    hub = connect_hub()
    session = hub.auth.sign_in(
        email=user_settings.user_email, password=user_settings.user_secret
    )
    hub.postgrest.auth(session.access_token)
    return hub


def push_instance():
    dobject_df = load("dobject")
    for id, v in dobject_df.index:
        push_dobject(id, v)


def push_dobject(id: str, v: str):
    instance = get_or_create_instance()
    dobject = load("dobject").loc[[(id, v)]].reset_index().to_dict("records")[0]
    if not jupynb_exists(dobject["jupynb_id"], dobject["jupynb_v"]):
        insert_jupynb(dobject["jupynb_id"], dobject["jupynb_v"], instance["id"])
    if not dobject_exists(id, v):
        insert_dobject(id, v)
    if not user_dobject_exists(id, v):
        insert_user_dobject(id, v)


def unpush_instance():
    instance = get_or_create_instance()
    dobject_df = load("dobject")
    jupynb_df = load("jupynb")
    for id, v in dobject_df.index:
        delete_user_dobject(id, v)
        delete_dobject(id, v)
    for id, v in jupynb_df.index:
        delete_jupynb(id, v)
    delete_user_instance(instance["id"])
    delete_instance()


# def get_user(user_id):
#     hub = get_hub_with_authentication()
#     data = (
#         hub.table("usermeta")
#         .select("*")
#         .eq("lnid", user_id)
#         .execute()
#     )
#     if len(data.data) > 0:
#         return data.data[0]
#     return None

# Instance


def get_or_create_instance():
    instance = get_instance()
    if instance is None:
        instance = insert_instance()
    if not user_instance_exists(instance["id"]):
        insert_user_instance(instance["id"])
    return instance


def get_instance():
    hub = get_hub_with_authentication()
    instance_settings = load_or_create_instance_settings()
    data = (
        hub.table("instance")
        .select("*")
        .eq("name", instance_settings.instance_name)
        .execute()
    )
    if len(data.data) > 0:
        return data.data[0]
    return None


def delete_instance():
    hub = get_hub_with_authentication()
    instance_settings = load_or_create_instance_settings()
    (
        hub.table("instance")
        .delete()
        .eq("name", instance_settings.instance_name)
        .execute()
    )


def insert_instance():
    hub = get_hub_with_authentication()
    instance_settings = load_or_create_instance_settings()
    data = (
        hub.table("instance")
        .insert(
            {
                "id": str(uuid.uuid4()),
                "name": instance_settings.instance_name,
                "owner_id": hub.auth.session().user.id.hex,
                "storage_dir": str(instance_settings.storage_dir),
                "dbconfig": instance_settings._dbconfig,
                "cache_dir": str(instance_settings.cache_dir),
                "sqlite_file": str(instance_settings._sqlite_file),
                "sqlite_file_local": str(instance_settings._sqlite_file_local),
                "db": instance_settings.db,
            }
        )
        .execute()
    )
    assert len(data.data) == 1
    return data.data[0]


# User instance


def user_instance_exists(instance_id):
    hub = get_hub_with_authentication()
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


def delete_user_instance(instance_id):
    hub = get_hub_with_authentication()
    (
        hub.table("user_instance")
        .delete()
        .eq("user_id", hub.auth.session().user.id.hex)
        .eq("instance_id", instance_id)
        .execute()
    )


def insert_user_instance(instance_id):
    hub = get_hub_with_authentication()
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


# Jupynb


def jupynb_exists(id, v):
    hub = get_hub_with_authentication()
    data = hub.table("jupynb").select("*").eq("id", id).eq("v", v).execute()
    if len(data.data) > 0:
        return True
    return False


def delete_jupynb(id, v):
    hub = get_hub_with_authentication()
    hub.table("jupynb").delete().eq("id", id).eq("v", v).execute()


def insert_jupynb(id, v, instance_id):
    hub = get_hub_with_authentication()
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


# Dobject


def dobject_exists(id, v):
    hub = get_hub_with_authentication()
    data = hub.table("dobject").select("*").eq("id", id).eq("v", v).execute()
    if len(data.data) > 0:
        return True
    return False


def delete_dobject(id, v):
    hub = get_hub_with_authentication()
    hub.table("dobject").delete().eq("id", id).eq("v", v).execute()


def insert_dobject(id, v):
    hub = get_hub_with_authentication()
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


def user_dobject_exists(id, v):
    hub = get_hub_with_authentication()
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


def delete_user_dobject(id, v):
    hub = get_hub_with_authentication()
    (
        hub.table("user_dobject")
        .delete()
        .eq("dobject_id", id)
        .eq("dobject_v", v)
        .eq("user_id", hub.auth.session().user.id.hex)
        .execute()
    )


def insert_user_dobject(id, v):
    hub = get_hub_with_authentication()
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
