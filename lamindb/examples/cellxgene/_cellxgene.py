from __future__ import annotations

from typing import TYPE_CHECKING, Collection, Literal, NamedTuple

if TYPE_CHECKING:
    from lamindb.base.types import FieldAttr
    from lamindb.models import Registry, Schema

CELLxGENEOrganisms = Literal[
    "human",
    "mouse",
    "zebra danio",
    "rhesus macaquedomestic pig",
    "chimpanzee",
    "white-tufted-ear marmoset",
    "sars-2",
]
FieldType = Literal["ontology_id", "name"]


def save_cellxgene_defaults() -> None:
    """Save default values of the CELLxGENE schema to the instance.

    Adds CELLxGENE specific (control) values that are not available in the ontologies:

    - "normal" Disease
    - "na" Ethnicity
    - "unknown" entries for DevelopmentalStage, Phenotype, and CellType
    - "tissue", "organoid", "primary cell culture", and "cell line" ULabels (tissue_type)
    - "cell", "nucleus", "na" ULabels (suspension_type)
    """
    import bionty as bt

    from lamindb.models import ULabel

    # "normal" in Disease
    normal = bt.Phenotype.from_source(
        ontology_id="PATO:0000461",
        source=bt.Source.get(name="pato", currently_used=True),
    )
    bt.Disease(
        uid=normal.uid,
        name=normal.name,
        ontology_id=normal.ontology_id,
        description=normal.description,
        source=normal.source,  # not sure
    ).save()

    # na, unknown
    for model, name in zip(
        [
            bt.Ethnicity,
            bt.Ethnicity,
            bt.DevelopmentalStage,
            bt.Phenotype,
            bt.CellType,
        ],
        ["na", "unknown", "unknown", "unknown", "unknown"],
    ):
        model(ontology_id=name, name=name, description="From CellxGene schema.").save()

    # tissue_type
    tissue_type = ULabel(
        name="TissueType",
        is_type=True,
        description='From CellxGene schema. Is "tissue", "organoid", "primary cell culture", or "cell line".',
    ).save()
    for name in ["tissue", "organoid", "primary cell culture", "cell line"]:
        ULabel(name=name, type=tissue_type, description="From CellxGene schema.").save()

    # suspension_type
    suspension_type = ULabel(
        name="SuspensionType",
        is_type=True,
        description='From CellxGene schema. This MUST be "cell", "nucleus", or "na".',
    ).save()
    for name in ["cell", "nucleus", "na"]:
        ULabel(
            name=name, type=suspension_type, description="From CellxGene schema."
        ).save()

    # organisms
    taxonomy_ids = [
        "NCBITaxon:9606",  # Homo sapiens (Human)
        "NCBITaxon:10090",  # Mus musculus (House mouse)
        "NCBITaxon:9544",  # Macaca mulatta (Rhesus monkey)
        "NCBITaxon:9825",  # Sus scrofa domesticus (Domestic pig)
        "NCBITaxon:9598",  # Pan troglodytes (Chimpanzee)
        "NCBITaxon:9483",  # Callithrix jacchus (White-tufted-ear marmoset)
        "NCBITaxon:7955",  # Danio rerio (Zebrafish)
    ]
    for ontology_id in taxonomy_ids:
        bt.Organism.from_source(
            ontology_id=ontology_id,
            source=bt.Source.get(name="ncbitaxon", currently_used=True),
        ).save()


def create_cellxgene_schema(
    *,
    field_types: FieldType | Collection[FieldType] = "ontology_id",
    spatial_library_id: str | None = None,
    organism: CELLxGENEOrganisms = "human",
) -> Schema:
    """Generates a :class:`~lamindb.Schema` for a specific CELLxGENE schema version.

    Args:
        field_types: One or several of 'ontology_id', 'name'.
        organism: The organism of the Schema.
        library_id: Identifier for the spatial library.
            Specifying this value enables curation against spatial requirements.
    """
    import bionty as bt

    from lamindb.models import Feature, Schema, ULabel

    class CategorySpec(NamedTuple):
        field: str | FieldAttr | list[Registry]
        default: str | None
        needs_organism: bool = False

    categoricals_to_spec: dict[str, CategorySpec] = {
        "assay": CategorySpec(bt.ExperimentalFactor.name, None, False),
        "assay_ontology_term_id": CategorySpec(
            bt.ExperimentalFactor.ontology_id, None, False
        ),
        "cell_type": CategorySpec(bt.CellType.name, "unknown", False),
        "cell_type_ontology_term_id": CategorySpec(
            bt.CellType.ontology_id, None, False
        ),
        "development_stage": CategorySpec(bt.DevelopmentalStage.name, "unknown", True),
        "development_stage_ontology_term_id": CategorySpec(
            bt.DevelopmentalStage.ontology_id, None, True
        ),
        "disease": CategorySpec(bt.Disease.name, "normal", False),
        "disease_ontology_term_id": CategorySpec(bt.Disease.ontology_id, None, False),
        "self_reported_ethnicity": CategorySpec(bt.Ethnicity.name, "unknown", False),
        "self_reported_ethnicity_ontology_term_id": CategorySpec(
            bt.Ethnicity.ontology_id, None, False
        ),
        "sex": CategorySpec(bt.Phenotype.name, "unknown", False),
        "sex_ontology_term_id": CategorySpec(bt.Phenotype.ontology_id, None, False),
        "suspension_type": CategorySpec(ULabel.name, "cell", False),
        "tissue": CategorySpec(bt.Tissue.name, None, False),
        "tissue_ontology_term_id": CategorySpec(
            [bt.Tissue.ontology_id, bt.CellType.ontology_id], None, False
        ),
        "tissue_type": CategorySpec(ULabel.name, "tissue", False),
        "organism": CategorySpec(bt.Organism.scientific_name, None, False),
        "organism_ontology_term_id": CategorySpec(bt.Organism.ontology_id, None, False),
        "donor_id": CategorySpec(str, "unknown", False),
    }

    def _get_source_cat_filters(
        field: str | FieldAttr | type[Registry], *, needs_organism: bool | None = None
    ) -> dict | None:
        """Some ontology are organism specific and their Features therefore need a `cat_filter`."""
        if isinstance(field, str) or not needs_organism:
            return None
        registry = field.field.model if hasattr(field, "field") else field
        entity = f"bionty.{registry.__name__}"
        filters = {"entity": entity, "currently_used": True}
        if needs_organism:
            filters["organism"] = organism
        return {"source": bt.Source.filter(**filters).one()}

    field_types_set = (
        {field_types} if isinstance(field_types, str) else set(field_types)
    )
    if field_types_set == {"ontology_id"}:
        categoricals = {
            k: v.field
            for k, v in categoricals_to_spec.items()
            if k.endswith("_ontology_term_id") or k == "donor_id"
        }
    elif field_types_set == {"name"}:
        categoricals = {
            k: v.field
            for k, v in categoricals_to_spec.items()
            if not k.endswith("_ontology_term_id") and k != "donor_id"
        }
    elif field_types_set == {"name", "ontology_id"}:
        categoricals = {k: v.field for k, v in categoricals_to_spec.items()}
    else:
        raise ValueError(
            f"Invalid field_types: {field_types}. Must contain 'ontology_id', 'name', or both."
        )

    organism_fields = {"organism", "organism_ontology_term_id"}
    obs_categoricals = {
        k: v for k, v in categoricals.items() if k not in organism_fields
    }

    var_schema = Schema(
        name="var of CELLxGENE",
        index=Feature(
            name="var_index",
            dtype=bt.Gene.ensembl_gene_id,
            cat_filters=_get_source_cat_filters(
                bt.Gene.ensembl_gene_id, needs_organism=True
            ),
        ).save(),
        itype=Feature,
        features=[Feature(name="feature_is_filtered", dtype=bool).save()],
        dtype="DataFrame",
        coerce=True,
    ).save()

    obs_features = []
    for field in obs_categoricals:
        if field == "var_index":
            continue
        dtype = obs_categoricals[field]
        needs_organism = categoricals_to_spec[field].needs_organism

        cat_filters: dict | list[dict] | None
        if isinstance(dtype, list):
            cat_filters = (
                [
                    _get_source_cat_filters(d, needs_organism=needs_organism)
                    for d in dtype
                ]
                if needs_organism
                else None
            )
        elif not isinstance(dtype, str):
            cat_filters = _get_source_cat_filters(dtype, needs_organism=needs_organism)
        else:
            cat_filters = None

        obs_features.append(
            Feature(  # type: ignore
                name=field,
                dtype=dtype,
                default_value=categoricals_to_spec[field].default,
                cat_filters=cat_filters,  # type: ignore
            ).save()
        )

    for name in ["is_primary_data", "suspension_type", "tissue_type"]:
        obs_features.append(Feature(name=name, dtype=ULabel.name).save())

    obs_schema = Schema(
        name=f"obs of CELLxGENE of {field_types}",
        features=obs_features,
        otype="DataFrame",
        minimal_set=True,
        coerce=True,
    ).save()

    slots = {"var": var_schema, "obs": obs_schema}

    uns_categoricals = {k: v for k, v in categoricals.items() if k in organism_fields}

    uns_features = [
        Feature(
            name=field,
            dtype=uns_categoricals[field],
            default_value=categoricals_to_spec[field].default,
        ).save()
        for field in uns_categoricals
    ]

    uns_schema = Schema(
        name="uns of CELLxGENE version",
        features=uns_features,
        otype="DataFrame",
        minimal_set=True,
        coerce=True,
    ).save()

    slots["uns"] = uns_schema

    # Add spatial validation if library_id is provided
    if spatial_library_id:
        scalefactors_schema = Schema(
            name=f"scalefactors of spatial {spatial_library_id}",
            features=[
                Feature(name="spot_diameter_fullres", dtype=float).save(),
                Feature(name="tissue_hires_scalef", dtype=float).save(),
            ],
        ).save()

        spatial_schema = Schema(
            name="CELLxGENE spatial metadata",
            features=[
                Feature(
                    name="is_single",
                    dtype=bool,
                    description="True if dataset represents single spatial unit (tissue section for Visium, array for Slide-seqV2)",
                ).save()
            ],
        ).save()

        slots["uns:spatial"] = spatial_schema
        slots[f"uns:spatial:{spatial_library_id}:scalefactors"] = scalefactors_schema

    # Spatial library ID must be in the name
    # Otherwise, we have lookup side effects where other existing Spatial Library IDs make it into the Schema
    schema_name = f"CELLxGENE AnnData of {', '.join(field_types) if isinstance(field_types, list) else field_types}"
    if spatial_library_id:
        schema_name += f" ({spatial_library_id})"

    full_cxg_schema = Schema(
        name=schema_name,
        otype="AnnData",
        minimal_set=True,
        coerce=True,
        slots=slots,
    ).save()

    return full_cxg_schema
