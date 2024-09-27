from django.contrib.postgres.aggregates import ArrayAgg, JSONBAgg
from django.db import connection
from django.db.models import F, Func, OuterRef, Q, Subquery
from django.db.models.fields.related import ForeignKey, ManyToManyField
from django.db.models.fields.reverse_related import ManyToManyRel, ManyToOneRel
from django.db.models.functions import JSONObject
from lnschema_core.models import Artifact, Collection

from .schema import dict_related_model_to_related_name


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


def get_non_id_value(item):
    return next((v for k, v in item.items() if k != "id"), "Unknown")


class DictAgg(Func):
    function = "JSONB_OBJECT_AGG"
    output_field = JSONBAgg()


def get_artifact_with_related(
    artifact: Artifact,
    fks: bool = True,
    m2m: bool = True,
    links: bool = True,
) -> dict:
    """Fetch an artifact with its related data.

    Args:
        artifact (Artifact | Collection): The artifact to fetch.
        m2m (bool, optional): Whether to include many-to-many relationships. Defaults to True.
        links (bool, optional): Whether to include link tables. Defaults to True.

    Returns:
        result: the fetched data.
    """
    from lamindb._can_validate import get_name_field

    fields = artifact._meta.fields
    direct_fields = []
    foreign_key_fields = []
    for f in fields:
        if f.is_relation:
            foreign_key_fields.append(f.name)
        else:
            direct_fields.append(f.name)

    model = artifact.__class__

    m2m_relations = (
        []
        if not m2m
        else [
            v
            for v in dict_related_model_to_related_name(model).values()
            if not v.startswith("_")
        ]
    )
    link_tables = (
        []
        if not links
        else list(dict_related_model_to_related_name(model, links=True).values())
    )

    # Clear previous queries
    connection.queries_log.clear()

    # Create subqueries for foreign key relationships
    fk_subqueries = {}
    if fks:
        for fk in foreign_key_fields:
            name_field = get_name_field(get_related_model(model, fk))
            fk_subqueries[f"fkfield_{fk}"] = JSONObject(
                id=F(f"{fk}__id"),
                display=F(f"{fk}__{name_field}"),
            )

    # Create subqueries for many-to-many relationships
    m2m_subqueries = {}
    if m2m:
        for name in m2m_relations:
            related_model = get_related_model(model, name)
            name_field = get_name_field(related_model)
            m2m_subqueries[f"m2mfield_{name}"] = ArrayAgg(
                JSONObject(id=F(f"{name}__id"), name=F(f"{name}__{name_field}")),
                filter=Q(**{f"{name}__isnull": False}),
                distinct=True,
            )

    # Create subqueries for link tables
    exclude_link_table_fields = {
        "run",
        "created_at",
        "created_by",
        "label_ref_is_name",
        "feature_ref_is_name",
    }
    link_subqueries = {}
    if links:
        for link in link_tables:
            link_model = getattr(model, link).rel.related_model
            fields = [
                f.name
                for f in link_model._meta.fields
                if f.name not in exclude_link_table_fields and f.name != "id"
            ]
            link_subqueries[f"linkfield_{link}"] = Subquery(
                link_model.objects.filter(artifact=OuterRef("pk"))
                .annotate(data=JSONObject(id=F("id"), **{f: F(f) for f in fields}))
                .values("artifact")
                .annotate(json_agg=ArrayAgg("data"))
                .values("json_agg")
            )

    # Combine all annotations
    annotations = {**m2m_subqueries, **fk_subqueries, **link_subqueries}

    # Execute the query
    artifact_meta = (
        model.objects.using(artifact._state.db)
        .filter(uid=artifact.uid)
        .annotate(**annotations)
        .values(
            *direct_fields,
            *m2m_subqueries.keys(),
            *fk_subqueries.keys(),
            *link_subqueries.keys(),
        )
        .first()
    )

    # Restructure the result
    if artifact_meta:
        related_data: dict = {"m2m": {}, "fk": {}, "link": {}}
        for k, v in artifact_meta.items():
            if k.startswith("m2mfield_"):
                related_data["m2m"][k[9:]] = v
            elif k.startswith("fkfield_"):
                related_data["fk"][k[8:]] = v
            elif k.startswith("linkfield_"):
                related_data["link"][k[10:]] = v

        # map m2m data to dictionaries
        m2m_data = {}
        for k, v in related_data["m2m"].items():
            if v:
                m2m_data[k] = {item["id"]: get_non_id_value(item) for item in v}
        related_data["m2m"] = m2m_data

        return {
            **{name: artifact_meta[name] for name in direct_fields},
            "related_data": related_data,
        }
    return None
