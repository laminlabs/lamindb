import pandas as pd

from .._logger import logger
from ..meta import annotate


class FeatureModel:
    def __init__(self, entity_model) -> None:
        self._entity_model = entity_model

    @property
    def entity(self):
        """Correspond to the feature entity table."""
        return self._entity_model.entity

    @property
    def id_type(self):
        """Type of id used for curation."""
        return self._entity_model._id_field

    @property
    def species(self):
        """Species."""
        return self._entity_model.species

    @property
    def df(self):
        """Reference table."""
        return self._entity_model.df

    def curate(self, df: pd.DataFrame):
        if self.id_type in df.columns:
            return self._entity_model.curate(df=df, column=self.id_type)
        else:
            logger.warning(f"{self.id_type} column not found, using index as features.")
            return self._entity_model.curate(df=df, column=None)

    def ingest(self, dobject_id, df_curated):
        """Ingest features."""
        annotate_feature = getattr(annotate, self.entity)
        mapped_df = self.df.loc[df_curated.index[df_curated["__curated__"]]].copy()
        mapped_dict = {}
        for i, row in mapped_df.iterrows():
            mapped_dict[i] = pd.concat([row, pd.Series([i], index=[self.id_type])])

        annotate_feature(
            dobject_id=dobject_id,
            values=mapped_dict,
            species=self.species,
        )
