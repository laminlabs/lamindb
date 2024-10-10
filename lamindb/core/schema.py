from __future__ import annotations

import lamindb_setup as ln_setup
from django.db.models import ManyToManyField
from lamindb_setup._connect_instance import (
    get_owner_name_from_identifier,
    load_instance_settings,
)
from lamindb_setup.core._settings_store import instance_settings_file
from lnschema_core.models import Feature, FeatureSet, LinkORM, Record


def get_schemas_modules(instance: str | None) -> set[str]:
    if instance is None:
        return ln_setup.settings.instance.schema
    owner, name = get_owner_name_from_identifier(instance)
    settings_file = instance_settings_file(name, owner)
    schema = load_instance_settings(settings_file).schema
    shared_schema_modules = ln_setup.settings.instance.schema.intersection(schema)
    shared_schema_modules.add("core")
    return shared_schema_modules


def dict_schema_name_to_model_name(
    registry: type[Record], instance: str | None = None
) -> dict[str, Record]:
    schema_modules = get_schemas_modules(instance)
    d: dict = {
        i.related_model.__get_name_with_schema__(): i.related_model
        for i in registry._meta.related_objects
        if i.related_name is not None
        and i.related_model.__get_schema_name__() in schema_modules
    }
    d.update(
        {
            i.related_model.__get_name_with_schema__(): i.related_model
            for i in registry._meta.many_to_many
            if i.name is not None
            and i.related_model.__get_schema_name__() in schema_modules
        }
    )
    return d


def dict_related_model_to_related_name(
    registry: type[Record], links: bool = False, instance: str | None = None
) -> dict[str, str]:
    def include(model: Record):
        return not links != issubclass(model, LinkORM)

    schema_modules = get_schemas_modules(instance)

    related_objects = registry._meta.related_objects + registry._meta.many_to_many
    d: dict = {
        record.related_model.__get_name_with_schema__(): (
            record.related_name
            if not isinstance(record, ManyToManyField)
            else record.name
        )
        for record in related_objects
        if (
            record.name is not None
            and include(record.related_model)
            and record.related_model.__get_schema_name__() in schema_modules
        )
    }
    return d


def get_related_name(features_type: type[Record]) -> str:
    candidates = [
        field.related_name
        for field in FeatureSet._meta.related_objects
        if field.related_model == features_type
    ]
    if not candidates:
        raise ValueError(
            f"Can't create feature sets from {features_type.__name__} because it's not"
            " related to it!\nYou need to create a link model between FeatureSet and"
            " your Record in your custom schema.\nTo do so, add a"
            " line:\nfeature_sets = models.ManyToMany(FeatureSet,"
            " related_name='mythings')\n"
        )
    return candidates[0]
