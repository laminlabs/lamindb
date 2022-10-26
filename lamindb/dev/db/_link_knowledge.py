from typing import Union

import pandas as pd
from lamin_logger import logger

from lamindb.dev.db._add import add
from lamindb.dev.db._link import link
from lamindb.dev.db._select import select
from lamindb.knowledge import CellMarker, Gene, Protein
from lamindb.schema import bionty


class LinkFeatureToKnowledgeTable:
    """Link features to a knowledge table."""

    def __init__(
        self,
        df: pd.DataFrame,
        knowledge_table: Union[Gene, Protein, CellMarker],
        featureset_name: str = None,
    ) -> None:
        self._model = knowledge_table
        self._id = knowledge_table._id_field
        self._featureset_name = featureset_name

        # Add and curate features against a knowledge table
        column = None
        if self._id in df.columns:
            column = self._id
        else:
            logger.warning(f"{self._id} column not found, using index as features.")
        df_curated = self._model.curate(df=df, column=column)

        # logging of curation
        n = df_curated["__curated__"].count()
        n_mapped = df_curated["__curated__"].sum()
        self._log = {
            "feature": self._id,
            "n_mapped": n_mapped,
            "percent_mapped": round(n_mapped / n * 100, 1),
            "unmapped": df_curated.index[~df_curated["__curated__"]],
        }
        self._df_curated = df_curated

    @property
    def log(self):
        """Logging."""
        return self._log

    def commit(self, dobject_id: str) -> None:
        """Commit features."""
        # mapped features will also contain fields in the reference table
        mapped_index = self._df_curated.index[self._df_curated["__curated__"]]
        mapped_df = self._model.df.loc[mapped_index].copy()
        mapped_dict = {}
        for i, row in mapped_df.iterrows():
            mapped_dict[i] = pd.concat([row, pd.Series([i], index=[self._id])])

        # unmapped features will only contain it's own field
        unmapped_dict = {}
        for um in self._df_curated.index.difference(mapped_index):
            unmapped_dict[um] = {self._id: um}

        species = select(bionty.species, common_name=self._model.species).one_or_none()
        if species is None:
            species = add(bionty.species(common_name=self._model.species))

        link.feature(
            dobject_id=dobject_id,
            features={**mapped_dict, **unmapped_dict},
            feature_entity=self._model.entity,
            species=self._model.species,
            featureset_name=self._featureset_name,
        )
