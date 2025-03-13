from pathlib import Path

import pandas as pd
import yaml  # type: ignore
from lamin_utils import logger


def _read_schema_versions(ontology_versions: Path) -> dict[str, pd.DataFrame]:
    data = yaml.safe_load(open(ontology_versions))
    schema_versions = data["schema-version"]

    def _schema_to_df(schema_data):
        return pd.DataFrame(
            [
                (entity, organism, ontology, version)
                for entity, details in schema_data.items()
                for ontology, values in details.items()
                for organism, version in values.items()
            ],
            columns=["entity", "source", "organism", "version"],
        ).set_index("entity")

    schema_versions_df = {
        version: _schema_to_df(details) for version, details in schema_versions.items()
    }

    return schema_versions_df


def _init_categoricals_additional_values(controls_were_created: bool | None = None):
    """Add additional values from CellxGene schema to the DB."""
    import bionty as bt

    import lamindb as ln

    # Note: if you add another control below, be mindful to change the if condition that
    # triggers whether creating these records is re-considered
    if controls_were_created is None:
        controls_were_created = (
            ln.ULabel.filter(name="SuspensionType", is_type=True).one_or_none()
            is not None
        )
    if not controls_were_created:
        logger.important("Creating control labels in the CellxGene schema.")
        bt.CellType(
            ontology_id="unknown", name="unknown", description="From CellxGene schema."
        ).save()
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
        bt.Ethnicity(
            ontology_id="na", name="na", description="From CellxGene schema."
        ).save()
        bt.Ethnicity(
            ontology_id="unknown", name="unknown", description="From CellxGene schema."
        ).save()
        bt.DevelopmentalStage(
            ontology_id="unknown", name="unknown", description="From CellxGene schema."
        ).save()
        bt.Phenotype(
            ontology_id="unknown", name="unknown", description="From CellxGene schema."
        ).save()

        # tissue_type
        tissue_type = ln.ULabel(
            name="TissueType",
            is_type=True,
            description='From CellxGene schema. Is "tissue", "organoid", or "cell culture".',
        ).save()
        ln.ULabel(
            name="tissue", type=tissue_type, description="From CellxGene schema."
        ).save()
        ln.ULabel(
            name="organoid", type=tissue_type, description="From CellxGene schema."
        ).save()
        ln.ULabel(
            name="cell culture", type=tissue_type, description="From CellxGene schema."
        ).save()

        # suspension_type
        suspension_type = ln.ULabel(
            name="SuspensionType",
            is_type=True,
            description='From CellxGene schema. This MUST be "cell", "nucleus", or "na".',
        ).save()
        ln.ULabel(
            name="cell", type=suspension_type, description="From CellxGene schema."
        ).save()
        ln.ULabel(
            name="nucleus", type=suspension_type, description="From CellxGene schema."
        ).save()
        ln.ULabel(name="na", type=suspension_type).save()
