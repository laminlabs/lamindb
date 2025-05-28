from __future__ import annotations

from functools import reduce
from typing import TYPE_CHECKING

from django.contrib.postgres.aggregates import ArrayAgg
from django.db import connection
from django.db.models import F, OuterRef, Q, Subquery
from django.db.models.fields.related import ForeignKey, ManyToManyField
from django.db.models.fields.reverse_related import ManyToManyRel, ManyToOneRel
from django.db.models.functions import JSONObject

from ._relations import dict_related_model_to_related_name, get_schema_modules
from .schema import Schema

if TYPE_CHECKING:
    from .artifact import Artifact
    from .sqlrecord import SQLRecord


def patch_many_to_many_descriptor() -> None:
    """Patches Django's `ManyToManyDescriptor.__get__` method to suggest better errors when saving relationships of an unsaved model.

    Before this patch: Cryptic errors are raised when relationships of an unsaved
    record are attempted to be modified.

    After this patch: Attempts to access M2M relationships on unsaved objects
    will raise ValueError, suggesting explicit .save() of the record to be modified
    before relationship creation.
    """
    from django.db.models.fields.related_descriptors import ManyToManyDescriptor

    original_get = ManyToManyDescriptor.__get__

    def patched_get(self, instance, cls=None):
        if instance is not None and instance.pk is None:
            raise ValueError(
                f"You are trying to access the many-to-many relationships of an unsaved {instance.__class__.__name__} object. Please save it first using '.save()'."
            )
        return original_get(self, instance, cls)

    ManyToManyDescriptor.__get__ = patched_get


def get_related_model(model, field_name):
    try:
        field = model._meta.get_field(field_name)

        if isinstance(field, (ForeignKey, ManyToManyField)):
            # Forward ForeignKey or ManyToManyField
            return field.remote_field.model
        elif isinstance(field, (ManyToOneRel, ManyToManyRel)):
            # Reverse ForeignKey or ManyToManyField
            return field.related_model
        else:
            return f"Unexpected field type: {type(field)}"
    except Exception as e:
        return f"Error: {str(e)}"


def get_artifact_with_related(
    artifact: SQLRecord,
    include_fk: bool = False,
    include_m2m: bool = False,
    include_feature_link: bool = False,
    include_schema: bool = False,
) -> dict:
    """Fetch an artifact with its related data."""
    from ._label_manager import EXCLUDE_LABELS
    from .can_curate import get_name_field

    model = artifact.__class__
    schema_modules = get_schema_modules(artifact._state.db)

    foreign_key_fields = [
        f.name
        for f in model._meta.fields
        if f.is_relation
        and f.related_model.__get_module_name__() in schema_modules
        and f.name != "branch"  # TODO: re-enable at some point
    ]

    m2m_relations = (
        []
        if not include_m2m
        else [
            v
            for v in dict_related_model_to_related_name(
                model, instance=artifact._state.db
            ).values()
            if not v.startswith("_") and v not in EXCLUDE_LABELS
        ]
    )
    link_tables = (
        []
        if not include_feature_link
        else list(
            dict_related_model_to_related_name(
                model, links=True, instance=artifact._state.db
            ).values()
        )
    )

    # Clear previous queries
    connection.queries_log.clear()

    annotations = {}

    if include_fk:
        for fk in foreign_key_fields:
            name_field = get_name_field(get_related_model(model, fk))
            annotations[f"fkfield_{fk}"] = JSONObject(
                id=F(f"{fk}__id"), name=F(f"{fk}__{name_field}")
            )

    for link in link_tables:
        link_model = getattr(model, link).rel.related_model
        if (
            not hasattr(link_model, "feature")
            or link_model.__name__ == "RecordArtifact"
        ):
            continue
        label_field = link.removeprefix("links_").replace("_", "")
        annotations[f"linkfield_{link}"] = Subquery(
            link_model.objects.filter(artifact=OuterRef("pk"))
            .annotate(
                data=JSONObject(
                    id=F("id"),
                    feature=F("feature"),
                    **{label_field: F(label_field)},
                )
            )
            .values("artifact")
            .annotate(json_agg=ArrayAgg("data"))
            .values("json_agg")
        )

    if include_schema:
        annotations["schemas"] = Subquery(
            model.feature_sets.through.objects.filter(artifact=OuterRef("pk"))
            .annotate(
                data=JSONObject(
                    id=F("id"),
                    slot=F("slot"),
                    schema=F("schema"),
                )
            )
            .values("artifact")
            .annotate(json_agg=ArrayAgg("data"))
            .values("json_agg")
        )

    artifact_meta = (
        model.objects.using(artifact._state.db)
        .filter(uid=artifact.uid)
        .annotate(**annotations)
        .values(*["id", "uid"], *annotations.keys())
        .first()
    )

    if not artifact_meta:
        return None

    related_data: dict = {"m2m": {}, "fk": {}, "link": {}, "schemas": {}}
    for k, v in artifact_meta.items():
        if k.startswith("fkfield_"):
            related_data["fk"][k[8:]] = v
        elif k.startswith("linkfield_"):
            related_data["link"][k[10:]] = v
        elif k == "schemas":
            if v:
                related_data["schemas"] = get_schema_m2m_relations(
                    artifact, {i["schema"]: i["slot"] for i in v}
                )

    if len(m2m_relations) == 0:
        m2m_any = False
    else:
        m2m_any_expr = reduce(
            lambda a, b: a | b,
            (Q(**{f"{m2m_name}__isnull": False}) for m2m_name in m2m_relations),
        )
        # this is needed to avoid querying all m2m relations even if they are all empty
        # this checks if non-empty m2m relations are present in the record
        m2m_any = (
            model.objects.using(artifact._state.db)
            .filter(uid=artifact.uid)
            .filter(m2m_any_expr)
            .exists()
        )
    if m2m_any:
        m2m_data = related_data["m2m"]
        for m2m_name in m2m_relations:
            related_model = get_related_model(model, m2m_name)
            name_field = get_name_field(related_model)
            m2m_records = (
                getattr(artifact, m2m_name).values_list("id", name_field).distinct()
            )
            for rec_id, rec_name in m2m_records:
                if m2m_name not in m2m_data:
                    m2m_data[m2m_name] = {}
                m2m_data[m2m_name][rec_id] = rec_name

    return {
        **{name: artifact_meta[name] for name in ["id", "uid"]},
        "related_data": related_data,
    }


def get_schema_m2m_relations(artifact: Artifact, slot_schema: dict, limit: int = 20):
    """Fetch all many-to-many relationships for given feature sets."""
    from .can_curate import get_name_field

    m2m_relations = [
        v
        for v in dict_related_model_to_related_name(Schema).values()
        if v is not None and not v.startswith("_") and v != "artifacts"
    ]

    annotations = {}
    related_names = {}
    for name in m2m_relations:
        related_model = get_related_model(Schema, name)
        if related_model is Schema:
            # this is for the `type` field
            continue
        name_field = get_name_field(related_model)

        # Get the correct field names for the through table
        if not hasattr(getattr(Schema, name), "through"):
            continue
        through_model = getattr(Schema, name).through

        # Subquery to get limited related records
        limited_related = Subquery(
            through_model.objects.filter(schema=OuterRef("pk")).values(
                related_model.__name__.lower()
            )[:limit]
        )

        annotations[f"m2mfield_{name}"] = ArrayAgg(
            JSONObject(id=F(f"{name}__id"), name=F(f"{name}__{name_field}")),
            filter=Q(
                **{
                    f"{name}__id__in": limited_related,
                }
            ),
            distinct=True,
        )
        related_names[name] = related_model.__get_name_with_module__()

    schema_m2m = (
        Schema.objects.using(artifact._state.db)
        .filter(id__in=slot_schema.keys())
        .annotate(**annotations)
        .values("id", *annotations.keys())
    )

    result = {}
    for fs in schema_m2m:
        slot = slot_schema.get(fs["id"])
        result[fs["id"]] = (
            slot,
            {
                related_names.get(k[9:]): [item["name"] for item in v]
                for k, v in fs.items()
                if k.startswith("m2mfield_") and v
            },
        )

    return result


patch_many_to_many_descriptor()
