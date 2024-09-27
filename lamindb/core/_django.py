from django.contrib.postgres.aggregates import ArrayAgg
from django.db import connection
from django.db.models import F, OuterRef, Q, Subquery
from django.db.models.fields.related import ForeignKey, ManyToManyField
from django.db.models.fields.reverse_related import ManyToManyRel, ManyToOneRel
from django.db.models.functions import JSONObject
from lnschema_core.models import Artifact, FeatureSet, Record

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


def get_artifact_with_related(
    artifact: Record,
    include_fk: bool = False,
    include_m2m: bool = False,
    include_feature_link: bool = False,
    include_featureset: bool = False,
) -> dict:
    """Fetch an artifact with its related data."""
    from lamindb._can_validate import get_name_field

    model = artifact.__class__
    foreign_key_fields = [f.name for f in model._meta.fields if f.is_relation]

    m2m_relations = (
        []
        if not include_m2m
        else [
            v
            for v in dict_related_model_to_related_name(model).values()
            if not v.startswith("_")
        ]
    )
    link_tables = (
        []
        if not include_feature_link
        else list(dict_related_model_to_related_name(model, links=True).values())
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

    for name in m2m_relations:
        related_model = get_related_model(model, name)
        name_field = get_name_field(related_model)
        annotations[f"m2mfield_{name}"] = ArrayAgg(
            JSONObject(id=F(f"{name}__id"), name=F(f"{name}__{name_field}")),
            filter=Q(**{f"{name}__isnull": False}),
            distinct=True,
        )

    for link in link_tables:
        link_model = getattr(model, link).rel.related_model
        if not hasattr(link_model, "feature"):
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

    if include_featureset:
        annotations["featuresets"] = Subquery(
            model.feature_sets.through.objects.filter(artifact=OuterRef("pk"))
            .annotate(
                data=JSONObject(
                    id=F("id"),
                    slot=F("slot"),
                    featureset=F("featureset"),
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

    related_data: dict = {"m2m": {}, "fk": {}, "link": {}, "featuresets": {}}
    for k, v in artifact_meta.items():
        if k.startswith("m2mfield_"):
            related_data["m2m"][k[9:]] = v
        elif k.startswith("fkfield_"):
            related_data["fk"][k[8:]] = v
        elif k.startswith("linkfield_"):
            related_data["link"][k[10:]] = v
        elif k == "featuresets":
            if v:
                related_data["featuresets"] = get_featureset_m2m_relations(
                    artifact, {i["featureset"]: i["slot"] for i in v}
                )

    related_data["m2m"] = {
        k: {item["id"]: item["name"] for item in v}
        for k, v in related_data["m2m"].items()
        if v
    }

    return {
        **{name: artifact_meta[name] for name in ["id", "uid"]},
        "related_data": related_data,
    }


def get_featureset_m2m_relations(
    artifact: Artifact, slot_featureset: dict, limit: int = 20
):
    """Fetch all many-to-many relationships for given feature sets."""
    from lamindb._can_validate import get_name_field

    m2m_relations = [
        v
        for v in dict_related_model_to_related_name(FeatureSet).values()
        if not v.startswith("_") and v != "artifacts"
    ]

    annotations = {}
    related_names = {}
    for name in m2m_relations:
        related_model = get_related_model(FeatureSet, name)
        name_field = get_name_field(related_model)

        # Get the correct field names for the through table
        through_model = getattr(FeatureSet, name).through
        related_field = (
            through_model.__name__.replace("FeatureSet", "").lower().replace("_", "")
        )

        # Subquery to get limited related records
        limited_related = Subquery(
            through_model.objects.filter(featureset=OuterRef("pk")).values(
                related_field
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
        related_names[name] = related_model.__get_name_with_schema__()

    featureset_m2m = (
        FeatureSet.objects.using(artifact._state.db)
        .filter(id__in=slot_featureset.keys())
        .annotate(**annotations)
        .values("id", *annotations.keys())
    )

    result = {}
    for fs in featureset_m2m:
        slot = slot_featureset.get(fs["id"])
        result[fs["id"]] = (
            slot,
            {
                related_names.get(k[9:]): [item["name"] for item in v]
                for k, v in fs.items()
                if k.startswith("m2mfield_") and v
            },
        )

    return result
