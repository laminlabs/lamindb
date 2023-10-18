from typing import Type
from uuid import UUID

import dj_database_url
from django.db import connections
from lamindb_setup._init_instance import InstanceSettings
from lamindb_setup._load_instance import get_owner_name_from_identifier
from lamindb_setup.dev._hub_core import load_instance
from lnschema_core import Registry

from lamindb._query_set import QuerySet


def add_db_connection(isettings: InstanceSettings, using: str):
    db_config = dj_database_url.config(
        default=isettings.db, conn_max_age=600, conn_health_checks=True
    )
    db_config["TIME_ZONE"] = "UTC"
    db_config["OPTIONS"] = {}
    db_config["AUTOCOMMIT"] = True

    connections.settings[using] = db_config


def filter(Registry: Type[Registry], using: str = None, **expressions) -> QuerySet:
    """See :meth:`~lamindb.dev.Registry.filter`."""
    if using is not None:
        owner, name = get_owner_name_from_identifier(using)
        load_result = load_instance(owner=owner, name=name)
        if isinstance(load_result, str):
            raise RuntimeError(
                f"Fail to load instance {using}, please check your permission!"
            )
        instance_result, storage_result = load_result
        isettings = InstanceSettings(
            owner=owner,
            name=name,
            storage_root=storage_result["root"],
            storage_region=storage_result["region"],
            db=instance_result["db"],
            schema=instance_result["schema_str"],
            id=UUID(instance_result["id"]),
        )
        add_db_connection(isettings, using)
    qs = QuerySet(model=Registry, using=using)
    if len(expressions) > 0:
        return qs.filter(**expressions)
    else:
        return qs
