from __future__ import annotations

from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from django.db.models import ManyToManyField
from lamindb_setup._connect_instance import (
    get_owner_name_from_identifier,
    load_instance_settings,
)
from lamindb_setup.core._settings_store import instance_settings_file

from lamindb.models.sqlrecord import IsLink

if TYPE_CHECKING:
    from lamindb.models.sqlrecord import Registry, SQLRecord


def get_schema_modules(instance: str | None) -> set[str]:
    if instance is None or instance == "default":
        schema_modules = set(ln_setup.settings.instance.modules)
        schema_modules.add("core")
        return schema_modules
    owner, name = get_owner_name_from_identifier(instance)
    settings_file = instance_settings_file(name, owner)
    if settings_file.exists():
        modules = set(load_instance_settings(settings_file).modules)
    else:
        cache_filepath = (
            ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
        )
        if cache_filepath.exists():
            modules = set(cache_filepath.read_text().split("\n")[1].split(","))
        else:
            raise ValueError(f"Instance {instance} not found")
    shared_schema_modules = set(ln_setup.settings.instance.modules).intersection(
        modules
    )
    shared_schema_modules.add("core")
    return shared_schema_modules


# this function here should likely be renamed
# it maps the __get_name_with_module__() onto the actual model
def dict_module_name_to_model_name(
    registry: Registry, instance: str | None = None
) -> dict[str, Registry]:
    schema_modules = get_schema_modules(instance)
    d: dict = {
        i.related_model.__get_name_with_module__(): i.related_model
        for i in registry._meta.related_objects
        if i.related_name is not None
        and i.related_model.__get_module_name__() in schema_modules
    }
    d.update(
        {
            i.related_model.__get_name_with_module__(): i.related_model
            for i in registry._meta.many_to_many
            if i.name is not None
            and i.related_model.__get_module_name__() in schema_modules
        }
    )
    return d


def dict_related_model_to_related_name(
    registry: type[SQLRecord], links: bool = False, instance: str | None = None
) -> dict[str, str]:
    def include(model: SQLRecord):
        return not links != issubclass(model, IsLink)

    schema_modules = get_schema_modules(instance)

    related_objects = registry._meta.related_objects + registry._meta.many_to_many
    d: dict = {
        record.related_model.__get_name_with_module__(): (
            record.related_name
            if not isinstance(record, ManyToManyField)
            else record.name
        )
        for record in related_objects
        if (
            record.name is not None
            and include(record.related_model)
            and record.related_model.__get_module_name__() in schema_modules
        )
    }
    return d


def get_related_name(features_type: type[SQLRecord]) -> str:
    from lamindb.models.schema import Schema

    candidates = [
        field.related_name
        for field in Schema._meta.related_objects
        if field.related_model == features_type
    ]
    if not candidates:
        raise ValueError(
            f"Can't create feature sets from {features_type.__name__} because it's not"
            " related to it!\nYou need to create a link model between Schema and"
            " your SQLRecord in your custom module.\nTo do so, add a"
            " line:\n_feature_sets = models.ManyToMany(Schema,"
            " related_name='mythings')\n"
        )
    return candidates[0]
