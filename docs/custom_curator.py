import bionty as bt
import pandas as pd

from lamindb.core import DataFrameCurator
from lamindb_setup.core.types import UPathStr
from lnschema_core import Record
from lnschema_core.types import FieldAttr

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

# Curate against these specified sources
FIXED_SOURCES = {
    "disease": bt.Source.filter(
        entity="bionty.Disease", name="mondo", version="2023-04-04"
    ).one(),
    "developmental_stage": bt.Source.filter(
        entity="bionty.DevelopmentalStage", name="hsapdv", version="2020-03-10"
    ).one(),
    "phenotype": bt.Source.filter(
        entity="bionty.Phenotype", name="hp", version="2023-06-17", organism="human"
    ).one(),
}


class EHRCurator(DataFrameCurator):
    """Custom curation flow for electronic health record data."""

    def __init__(
        self,
        data: pd.DataFrame | UPathStr,
        categoricals: dict[str, FieldAttr] = DEFAULT_CATEGORICALS,
        *,
        defaults: dict[str, str] = None,
        sources: dict[str, Record] = FIXED_SOURCES,
        organism="human",
    ):
        self.data = data
        self.organism = organism

        if defaults:
            for col, default in defaults.items():
                if col not in data.columns:
                    data[col] = default
                else:
                    data[col].fillna(default, inplace=True)

        super().__init__(
            df=data, categoricals=categoricals, sources=sources, organism=organism
        )

    def validate(self, organism: str | None = None) -> bool:
        """Further custom validation."""
        # --- Custom validation logic goes here --- #
        return DataFrameCurator.validate(self, organism)
