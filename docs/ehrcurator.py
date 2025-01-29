import bionty as bt
import pandas as pd
from lamindb.core import DataFrameCatCurator, logger
from lamindb.core.types import UPathStr

__version__ = "0.1.0"


class EHRCurator(DataFrameCatCurator):
    """Custom curation flow for electronic health record data."""

    def __init__(self, data: pd.DataFrame | UPathStr):
        # Curate these columns against the specified fields
        DEFAULT_CATEGORICALS = {
            "disease": bt.Disease.name,
            "phenotype": bt.Phenotype.name,
            "developmental_stage": bt.DevelopmentalStage.name,
        }

        # If columns or values are missing, we substitute with these defaults
        DEFAULT_VALUES = {
            "disease": "normal",
            "development_stage": "unknown",
            "phenotype": "unknown",
        }

        # Validate values onto the following ontology versions
        DEFAULT_SOURCES = {
            "disease": bt.Source.get(
                entity="bionty.Disease", name="mondo", version="2023-04-04"
            ),
            "developmental_stage": bt.Source.get(
                entity="bionty.DevelopmentalStage", name="hsapdv", version="2020-03-10"
            ),
            "phenotype": bt.Source.get(
                entity="bionty.Phenotype",
                name="hp",
                version="2023-06-17",
                organism="human",
            ),
        }

        self.data = data

        for col, default in DEFAULT_VALUES.items():
            if col not in self.data.columns:
                self.data[col] = default
            else:
                self.data[col] = self.data[col].fillna(default)

        super().__init__(
            df=self.data,
            categoricals=DEFAULT_CATEGORICALS,
            sources=DEFAULT_SOURCES,
            organism="human",
        )

    def validate(self, organism: str | None = None) -> bool:
        """Validates the internal EHR standard."""
        missing_columns = {"disease", "phenotype", "developmental_stage", "age"} - set(
            self.data.columns
        )
        if missing_columns:
            logger.error(
                f"Columns {', '.join(map(repr, missing_columns))} are missing but required."
            )
            return False

        return DataFrameCatCurator.validate(self, organism)
