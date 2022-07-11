from lamindb._setup import (
    load_or_create_instance_settings,
    load_or_create_user_settings,
)

from .._setup._hub import connect_hub
from ._load import load


def push(
    dobject_id=None,
):
    """Push an instance or a single dobject if specified."""
    user_settings = load_or_create_user_settings()
    instance_settings = load_or_create_instance_settings()

    hub = connect_hub()
    session = hub.auth.sign_in(
        email=user_settings.user_email, password=user_settings.user_secret
    )
    hub.postgrest.auth(session.access_token)

    #

    instance = {
        "id": instance_settings.instance_name,
        "name": instance_settings.instance_name,
    }

    instance_exists = (
        hub.table("instance")
        .select("*", count="exact")
        .eq("id", instance["id"])
        .execute()
        .count
        > 0  # noqa
    )

    if not instance_exists:
        data = hub.table("instance").insert(instance).execute()
        assert len(data.data) > 0

    user_instance = {
        "instance_id": instance_settings.instance_name,
        "user_id": session.user.id.hex,
    }

    user_instance_exists = (
        hub.table("user_instance")
        .select("*", count="exact")
        .eq("instance_id", user_instance["instance_id"])
        .eq("user_id", user_instance["user_id"])
        .execute()
        .count
        > 0  # noqa
    )

    if not user_instance_exists:
        data = hub.table("user_instance").insert(user_instance).execute()
        assert len(data.data) > 0

    #

    jupynb_df = load("jupynb").reset_index()
    dobject_df = load("dobject").reset_index()

    jupynb_df["time_created"] = jupynb_df["time_created"].dt.strftime(
        "%Y-%m-%d %H:%M:%S"
    )
    jupynb_df["time_updated"] = jupynb_df["time_updated"].dt.strftime(
        "%Y-%m-%d %H:%M:%S"
    )
    jupynb_df = jupynb_df.drop("user_id", axis=1)

    dobject_df["time_created"] = dobject_df["time_created"].dt.strftime(
        "%Y-%m-%d %H:%M:%S"
    )
    dobject_df["time_updated"] = dobject_df["time_updated"].dt.strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    if dobject_id:
        push_object(instance["id"], dobject_id, dobject_df, jupynb_df, hub)
    else:
        for dobject_id in dobject_df["id"].values:
            push_object(instance["id"], dobject_id, dobject_df, jupynb_df, hub)


def push_object(instance_id, dobject_id, dobject_df, jupynb_df, hub):
    """Push a dobject into the hub."""
    dobject_exists = (
        hub.table("dobject")
        .select("*", count="exact")
        .eq("id", dobject_id)
        .execute()
        .count
        > 0  # noqa
    )

    if not dobject_exists:
        dobject = dobject_df[dobject_df["id"] == dobject_id].to_dict("records")[0]
        push_jupynb(instance_id, dobject["jupynb_id"], jupynb_df, hub)
        data = hub.table("dobject").insert(dobject).execute()
        assert len(data.data) > 0


def push_jupynb(instance_id, jupynb_id, jupynb_df, hub):
    """Push a jupynb into the hub."""
    jupynb_exists = (
        hub.table("jupynb")
        .select("*", count="exact")
        .eq("id", jupynb_id)
        .execute()
        .count
        > 0  # noqa
    )

    if not jupynb_exists:
        jupynb = jupynb_df[jupynb_df["id"] == jupynb_id].to_dict("records")[0]
        data = (
            hub.table("jupynb").insert({**jupynb, "instance_id": instance_id}).execute()
        )
        assert len(data.data) > 0


def unpush(
    dobject_id=None,
):
    """Unpush an instance or a single dobject if specified."""
    user_settings = load_or_create_user_settings()
    instance_settings = load_or_create_instance_settings()

    hub = connect_hub()
    session = hub.auth.sign_in(
        email=user_settings.user_email, password=user_settings.user_secret
    )
    hub.postgrest.auth(session.access_token)

    jupynb_df = load("jupynb").reset_index()
    dobject_df = load("dobject").reset_index()

    if dobject_id:
        unpush_object(dobject_id, hub)
    else:
        for dobject_id in dobject_df["id"].values:
            unpush_object(dobject_id, hub)
        for jupynb_id in jupynb_df["id"].values:
            unpush_jupynb(jupynb_id, hub)
        data = (
            hub.table("user_instance")
            .delete()
            .eq("instance_id", instance_settings.instance_name)
            .eq("user_id", session.user.id.hex)
            .execute()
        )
        assert len(data.data) > 0
        data = (
            hub.table("instance")
            .delete()
            .eq("id", instance_settings.instance_name)
            .execute()
        )
        assert len(data.data) > 0


def unpush_object(dobject_id, hub):
    """Reverse push of a dobject into the hub."""
    data = hub.table("dobject").delete().eq("id", dobject_id).execute()
    assert len(data.data) > 0


def unpush_jupynb(jupynb_id, hub):
    """Reverse push of a dobject into the hub."""
    data = hub.table("jupynb").delete().eq("id", jupynb_id).execute()
    assert len(data.data) > 0
