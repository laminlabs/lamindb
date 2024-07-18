from __future__ import annotations

from django.db.models import ManyToManyField
from lnschema_core.models import Feature, FeatureSet, LinkORM, Record


def dict_schema_name_to_model_name(orm: type[Record]) -> dict[str, Record]:
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


def dict_related_model_to_related_name(
    orm: type[Record], links: bool = False
) -> dict[str, str]:
    def include(model: Record):
        return not links != issubclass(model, LinkORM)

    related_objects = orm._meta.related_objects + orm._meta.many_to_many
    d: dict = {
        record.related_model.__get_name_with_schema__(): (
            record.related_name
            if not isinstance(record, ManyToManyField)
            else record.name
        )
        for record in related_objects
        if (record.name is not None and include(record.related_model))
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
