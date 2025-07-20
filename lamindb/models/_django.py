from __future__ import annotations

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

    Before this patch: Cryptic errors are raised when relationships of an unsaved record are attempted to be modified.

    After this patch: Attempts to access M2M relationships on unsaved objects will raise ValueError,
    suggesting explicit .save() of the record to be modified before relationship creation.
    """
    from django.db.models.fields.related_descriptors import ManyToManyDescriptor

    original_get = ManyToManyDescriptor.__get__

    def patched_get(self, instance, cls=None):
        if instance is not None and instance.pk is None:
            raise ValueError(
                f"You are trying to access the many-to-many relationships of an unsaved {instance.__class__.__name__} object. Please save it first using '.save()'."
            )

        manager = original_get(self, instance, cls)
        if manager is None or not hasattr(manager, "add"):
            return manager

        original_manager_add = manager.add

        def patched_manager_add(*objs, **kwargs):
            try:
                return original_manager_add(*objs, **kwargs)
            except ValueError as e:
                if "Cannot add" in str(e) and "database" in str(e):
                    source_db = manager.instance._state.db

                    raise ValueError(
                        f"Cannot label a record from instance '{source_db}'. "
                        f"Please save the record first to your instance using '.save()'."
                    ) from None
                raise

        manager.add = patched_manager_add
        return manager

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
        if f.is_relation and f.related_model.__get_module_name__() in schema_modules
    ]

    # Create the map that the conversion function will need.
    # It maps the target model class to the m2m field name, e.g.,
    # {<class 'Ulabel'>: 'ulabels', <class 'CellType'>: 'cell_types'}
    m2m_model_to_field_map = {}
    if include_m2m:
        full_map = dict_related_model_to_related_name(
            model, instance=artifact._state.db
        )
        m2m_model_to_field_map = {
            model_cls: field_name
            for model_cls, field_name in full_map.items()
            if not field_name.startswith("_") and field_name not in EXCLUDE_LABELS
        }
    list(m2m_model_to_field_map.values())
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
            if fk == "run":
                annotations[f"fkfield_{fk}"] = JSONObject(
                    id=F(f"{fk}__id"),
                    name=F(f"{fk}__{name_field}"),
                    transform_key=F(f"{fk}__transform__key"),
                )
            else:
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
        related_model = link_model._meta.get_field(label_field).related_model
        name_field = get_name_field(related_model)
        label_field_name = f"{label_field}__{name_field}"
        annotations[f"linkfield_{link}"] = Subquery(
            link_model.objects.filter(artifact=OuterRef("pk"))
            .annotate(
                data=JSONObject(
                    id=F("id"),
                    feature=F("feature"),
                    **{label_field: F(label_field)},
                    **{label_field + "_display": F(label_field_name)},
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
        if k.startswith("fkfield_") and v is not None:
            related_data["fk"][k[8:]] = v
        elif k.startswith("linkfield_") and v is not None:
            related_data["link"][k[10:]] = v
        elif k == "schemas":
            if v:
                related_data["schemas"] = get_schema_m2m_relations(
                    artifact, {i["schema"]: i["slot"] for i in v}
                )

    def convert_link_data_to_m2m(
        link_data: dict,
        model,  # The main artifact model class is still needed for introspection
        m2m_model_map: dict,  # The pre-computed map from Step 1
    ) -> dict:
        """Converts link data to M2M-style data using a pre-computed model-to-field-name map."""
        m2m_data = {}
        for link_name, records in link_data.items():
            if not records:
                continue
            link_model = getattr(model, link_name).rel.related_model
            id_field_name = link_name.removeprefix("links_").replace("_", "")
            final_target_model = link_model._meta.get_field(id_field_name).related_model
            m2m_field_name = m2m_model_map.get(
                final_target_model.__get_name_with_module__()
            )
            display_field_name = f"{id_field_name}_display"
            m2m_data[m2m_field_name] = {
                record[id_field_name]: record[display_field_name] for record in records
            }
        return m2m_data

    related_data["m2m"] = convert_link_data_to_m2m(
        related_data["link"], model=model, m2m_model_map=m2m_model_to_field_map
    )
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
