import pandas as pd
from lamin_utils import logger
from lamindb_setup.core.upath import UPath

from lamindb.base.types import FieldAttr
from lamindb.models import Record, ULabel
from lamindb.models._from_values import _format_values

RESERVED_NAMES = {
    "ethnicity",
    "ethnicity_ontology_term_id",
    "X_normalization",
    "default_field",
    "layer_descriptions",
    "tags",
    "versions",
    "contributors",
    "preprint_doi",
    "project_description",
    "project_links",
    "project_name",
    "publication_doi",
}


def _get_cxg_categoricals() -> dict[str, FieldAttr]:
    import bionty as bt

    return {
        "assay": bt.ExperimentalFactor.name,
        "assay_ontology_term_id": bt.ExperimentalFactor.ontology_id,
        "cell_type": bt.CellType.name,
        "cell_type_ontology_term_id": bt.CellType.ontology_id,
        "development_stage": bt.DevelopmentalStage.name,
        "development_stage_ontology_term_id": bt.DevelopmentalStage.ontology_id,
        "disease": bt.Disease.name,
        "disease_ontology_term_id": bt.Disease.ontology_id,
        # "donor_id": "str",  via pandera
        "self_reported_ethnicity": bt.Ethnicity.name,
        "self_reported_ethnicity_ontology_term_id": bt.Ethnicity.ontology_id,
        "sex": bt.Phenotype.name,
        "sex_ontology_term_id": bt.Phenotype.ontology_id,
        "suspension_type": ULabel.name,
        "tissue": bt.Tissue.name,
        "tissue_ontology_term_id": bt.Tissue.ontology_id,
        "tissue_type": ULabel.name,
        "organism": bt.Organism.name,
        "organism_ontology_term_id": bt.Organism.ontology_id,
    }


def _create_sources(
    categoricals: dict[str, FieldAttr], schema_version: str, organism: str
) -> dict[str, Record]:
    """Creates a sources dictionary that can be passed to AnnDataCatManager."""
    import bionty as bt

    def _fetch_bionty_source(entity: str, organism: str) -> Record | None:  # type: ignore
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
                logger.error(
                    f"Could not find source: {entity}\n"
                    "    → consider running `bionty.core.sync_all_sources_to_latest()` and re-connect to your instance"
                )
            return source

    sources_df = pd.read_csv(UPath(__file__).parent / "schema_versions.csv")
    sources_df = sources_df[sources_df.schema_version == schema_version]
    if sources_df.empty:
        raise ValueError(
            f"Invalid schema_version: {schema_version}\n"
            f"Valid versions are: {_format_values(sources_df.schema_version.unique())}"
        )

    key_to_source: dict[str, bt.Source] = {}
    for key, field in categoricals.items():
        if field.field.model.__get_module_name__() == "bionty":
            entity = field.field.model.__name__
            key_to_source[key] = _fetch_bionty_source(entity, organism)
    key_to_source["var_index"] = _fetch_bionty_source("Gene", organism)

    return key_to_source


def _init_categoricals_additional_values() -> None:
    """Add additional values from CellxGene schema to the DB."""
    import bionty as bt

    # Note: if you add another control below, be mindful to change the if condition that
    # triggers whether creating these records is re-considered
    controls_were_created = (
        ULabel.filter(name="SuspensionType", is_type=True).one_or_none() is not None
    )
    if not controls_were_created:
        logger.important("Creating control labels in the CellxGene schema.")

        # "normal" in Disease
        normal = bt.Phenotype.from_source(
            ontology_id="PATO:0000461",
            source=bt.Source.get(name="pato", version="2024-03-28"),
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
            model(
                ontology_id=name, name=name, description="From CellxGene schema."
            ).save()

        # tissue_type
        tissue_type = ULabel(
            name="TissueType",
            is_type=True,
            description='From CellxGene schema. Is "tissue", "organoid", or "cell culture".',
        ).save()
        for name in ["tissue", "organoid", "cell culture"]:
            ULabel(
                name=name, type=tissue_type, description="From CellxGene schema."
            ).save()

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
