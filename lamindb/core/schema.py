from typing import Type

from lnschema_core.models import FeatureSet, LinkORM, Registry


def dict_schema_name_to_model_name(orm: Type[Registry]) -> dict[str, Registry]:
    d: dict = {
        i.related_model.__get_name_with_schema__(): i.related_model
        for i in orm._meta.related_objects
        if i.related_name is not None
    }
    d.update(
        {
            i.related_model.__get_name_with_schema__(): i.related_model
            for i in orm._meta.many_to_many
            if i.name is not None
        }
    )
    return d


def dict_related_model_to_related_name(orm: Type[Registry]) -> dict[str, str]:
    d: dict = {
        i.related_model.__get_name_with_schema__(): i.related_name
        for i in orm._meta.related_objects
        if i.related_name is not None
    }
    d.update(
        {
            i.related_model.__get_name_with_schema__(): i.name
            for i in orm._meta.many_to_many
            if (i.name is not None and not issubclass(i.related_model, LinkORM))
        }
    )

    return d


def get_related_name(features_type: Registry) -> str:
    candidates = [
        field.related_name
        for field in FeatureSet._meta.related_objects
        if field.related_model == features_type
    ]
    if not candidates:
        raise ValueError(
            f"Can't create feature sets from {features_type.__name__} because it's not"
            " related to it!\nYou need to create a link model between FeatureSet and"
            " your Registry in your custom schema.\nTo do so, add a"
            " line:\nfeature_sets = models.ManyToMany(FeatureSet,"
            " related_name='mythings')\n"
        )
    return candidates[0]
