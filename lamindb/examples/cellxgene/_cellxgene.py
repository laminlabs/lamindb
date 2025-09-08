from __future__ import annotations

from typing import TYPE_CHECKING, Collection, Literal, NamedTuple

import pandas as pd
from lamindb_setup.core import deprecated
from lamindb_setup.core.upath import UPath
from packaging import version

from lamindb.models._from_values import _format_values

if TYPE_CHECKING:
    from lamindb.base.types import FieldAttr
    from lamindb.models import Schema, SQLRecord

CELLxGENESchemaVersions = Literal["4.0.0", "5.0.0", "5.1.0", "5.2.0", "5.3.0", "6.0.0"]
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


@deprecated(new_name="save_cellxgene_defaults")
def save_cxg_defaults() -> None:
    return save_cellxgene_defaults()


def save_cellxgene_defaults() -> None:
    """Save default values of the CELLxGENE schema to the instance.

    Adds CELLxGENE specific (control) values that are not available in the ontologies:

    - "normal" Disease
    - "na" Ethnicity
    - "unknown" entries for DevelopmentalStage, Phenotype, and CellType
    - "tissue", "organoid", and "cell culture" ULabels (tissue_type)
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
        description='From CellxGene schema. Is "tissue", "organoid", or "cell culture".',
    ).save()
    for name in ["tissue", "organoid", "cell culture"]:
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


def _create_cellxgene_sources(
    categoricals: dict[str, FieldAttr], schema_version: str, organism: str
) -> dict[str, SQLRecord]:
    """Create a source dictionary of CELLxGENE categoricals to Source."""
    import bionty as bt

    def _fetch_bionty_source(entity: str, organism: str) -> SQLRecord | None:  # type: ignore
        """Fetch the Bionty source of the pinned ontology."""
        entity_sources = sources_df.loc[(sources_df.entity == entity)].copy()
        if not entity_sources.empty:
            if len(entity_sources) == 1:
                row = entity_sources.iloc[0]  # for sources with organism "all"
            else:
                row = entity_sources[entity_sources.organism == organism].iloc[0]
            source = bt.Source.filter(
                organism=row.organism,
                entity=f"bionty.{entity}",
                name=row.source,
                version=row.version,
            ).one_or_none()
            if source is None:
                source = getattr(bt, entity).add_source(
                    source=row.source,
                    version=row.version,
                    organism=row.organism,
                )
            return source

    sources_df = pd.read_csv(UPath(__file__).parent / "cellxgene_schema_versions.csv")
    sources_df = sources_df[sources_df.schema_version == schema_version]
    if sources_df.empty:
        raise ValueError(
            f"Invalid schema_version: {schema_version}\n"
            f"Valid versions are: {_format_values(sources_df.schema_version.unique())}"
        )

    key_to_source: dict[str, bt.Source] = {}
    for key, field in categoricals.items():
        if hasattr(field, "field"):
            if field.field.model.__get_module_name__() == "bionty":
                entity = field.field.model.__name__
                key_to_source[key] = _fetch_bionty_source(entity, organism)
        else:
            key_to_source[key] = field
    key_to_source["var_index"] = _fetch_bionty_source("Gene", organism)

    return key_to_source


@deprecated(new_name="create_cellxgene_schema")
def get_cxg_schema(
    schema_version: CELLxGENESchemaVersions,
    *,
    field_types: FieldType | Collection[FieldType] = "ontology_id",
    organism: CELLxGENEOrganisms = "human",
    spatial_library_id: str | None = None,
) -> Schema:
    return create_cellxgene_schema(
        schema_version,
        field_types=field_types,
        organism=organism,
        spatial_library_id=spatial_library_id,
    )


def create_cellxgene_schema(
    schema_version: CELLxGENESchemaVersions,
    *,
    field_types: FieldType | Collection[FieldType] = "ontology_id",
    organism: CELLxGENEOrganisms = "human",
    spatial_library_id: str | None = None,
) -> Schema:
    """Generates a :class:`~lamindb.Schema` for a specific CELLxGENE schema version.

    Args:
        schema_version: The CELLxGENE Schema version.
        field_types: One or several of 'ontology_id', 'name'.
        organism: The organism of the Schema.
        library_id: Identifier for the spatial library.
            Specifying this value enables curation against spatial requirements.
    """
    import bionty as bt

    from lamindb.models import Feature, Schema, ULabel

    class CategorySpec(NamedTuple):
        field: str | FieldAttr
        default: str | None

    categoricals_to_spec: dict[str, CategorySpec] = {
        "assay": CategorySpec(bt.ExperimentalFactor.name, None),
        "assay_ontology_term_id": CategorySpec(bt.ExperimentalFactor.ontology_id, None),
        "cell_type": CategorySpec(bt.CellType.name, "unknown"),
        "cell_type_ontology_term_id": CategorySpec(bt.CellType.ontology_id, None),
        "development_stage": CategorySpec(bt.DevelopmentalStage.name, "unknown"),
        "development_stage_ontology_term_id": CategorySpec(
            bt.DevelopmentalStage.ontology_id, None
        ),
        "disease": CategorySpec(bt.Disease.name, "normal"),
        "disease_ontology_term_id": CategorySpec(bt.Disease.ontology_id, None),
        "self_reported_ethnicity": CategorySpec(bt.Ethnicity.name, "unknown"),
        "self_reported_ethnicity_ontology_term_id": CategorySpec(
            bt.Ethnicity.ontology_id, None
        ),
        "sex": CategorySpec(bt.Phenotype.name, "unknown"),
        "sex_ontology_term_id": CategorySpec(bt.Phenotype.ontology_id, None),
        "suspension_type": CategorySpec(ULabel.name, "cell"),
        "tissue": CategorySpec(bt.Tissue.name, None),
        "tissue_ontology_term_id": CategorySpec(bt.Tissue.ontology_id, None),
        "tissue_type": CategorySpec(ULabel.name, "tissue"),
        "organism": CategorySpec(bt.Organism.scientific_name, None),
        "organism_ontology_term_id": CategorySpec(bt.Organism.ontology_id, None),
        "donor_id": CategorySpec(str, "unknown"),
    }

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

    is_version_6_or_later = version.parse(schema_version) >= version.parse("6.0.0")

    organism_fields = {"organism", "organism_ontology_term_id"}
    if is_version_6_or_later:
        obs_categoricals = {
            k: v for k, v in categoricals.items() if k not in organism_fields
        }
    else:
        obs_categoricals = categoricals

    sources = _create_cellxgene_sources(
        categoricals=categoricals,
        schema_version=schema_version,
        organism=organism,
    )

    var_schema = Schema(
        name=f"var of CELLxGENE version {schema_version}",
        index=Feature(
            name="var_index",
            dtype=bt.Gene.ensembl_gene_id,
            cat_filters={"source": sources["var_index"]},
        ).save(),
        itype=Feature,
        features=[Feature(name="feature_is_filtered", dtype=bool).save()],
        dtype="DataFrame",
        coerce_dtype=True,
    ).save()

    obs_features = [
        Feature(
            name=field,
            dtype=obs_categoricals[field],
            cat_filters={"source": source},
            default_value=categoricals_to_spec[field].default,
        ).save()
        for field, source in sources.items()
        if field != "var_index" and field in obs_categoricals
    ]
    for name in ["is_primary_data", "suspension_type", "tissue_type"]:
        obs_features.append(Feature(name=name, dtype=ULabel.name).save())

    obs_schema = Schema(
        name=f"obs of CELLxGENE version {schema_version} for {organism} of {field_types}",
        features=obs_features,
        otype="DataFrame",
        minimal_set=True,
        coerce_dtype=True,
    ).save()

    slots = {"var": var_schema, "obs": obs_schema}

    if is_version_6_or_later:
        uns_categoricals = {
            k: v for k, v in categoricals.items() if k in organism_fields
        }

        uns_features = [
            Feature(
                name=field,
                dtype=uns_categoricals[field],
                cat_filters={"source": sources[field]},
                default_value=categoricals_to_spec[field].default,
            ).save()
            for field in uns_categoricals
        ]

        uns_schema = Schema(
            name=f"uns of CELLxGENE version {schema_version}",
            features=uns_features,
            otype="DataFrame",
            minimal_set=True,
            coerce_dtype=True,
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
            slots[f"uns:spatial:{spatial_library_id}:scalefactors"] = (
                scalefactors_schema
            )

    full_cxg_schema = Schema(
        name=f"AnnData of CELLxGENE version {schema_version} for {organism} of {', '.join(field_types) if isinstance(field_types, list) else field_types}",
        otype="AnnData",
        minimal_set=True,
        coerce_dtype=True,
        slots=slots,
    ).save()

    return full_cxg_schema
