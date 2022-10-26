from typing import Union

import pandas as pd
from bionty import CellMarker, Gene, Protein
from lamin_logger import logger

from lamindb.dev.db._link import link


def add_features_from_knowledge_table(
    df: pd.DataFrame,
    knowledge_table: Union[Gene, Protein, CellMarker],
    featureset_name: str = None,
):
    """Curate a DataFrame with a feature model."""
    fm = LinkFeatureToKnowledgeTable(knowledge_table, featureset_name=featureset_name)
    df_curated, log = fm.curate(df)

    return {"model": fm, "df_curated": df_curated, "log": log}


class LinkFeatureToKnowledgeTable:
    def __init__(
        self,
        knowledge_table: Union[Gene, Protein, CellMarker],
        featureset_name: str = None,
    ) -> None:
        self._model = knowledge_table
        self._id = knowledge_table._id_field
        self._featureset_name = featureset_name

    def curate(self, df: pd.DataFrame):
        column = None
        if self._id in df.columns:
            column = self._id
        else:
            logger.warning(f"{self._id} column not found, using index as features.")
        df_curated = self._model.curate(df=df, column=column)

        # logging of curation
        n = df_curated["__curated__"].count()
        n_mapped = df_curated["__curated__"].sum()
        log = {
            "feature": self._id,
            "n_mapped": n_mapped,
            "percent_mapped": round(n_mapped / n * 100, 1),
            "unmapped": df_curated.index[~df_curated["__curated__"]],
        }

        return df_curated, log

    def ingest(self, dobject_id: str, df_curated: pd.DataFrame):
        """Ingest features."""
        # mapped features will also contain fields in the reference table
        mapped_index = df_curated.index[df_curated["__curated__"]]
        mapped_df = self._model.df.loc[mapped_index].copy()
        mapped_dict = {}
        for i, row in mapped_df.iterrows():
            mapped_dict[i] = pd.concat([row, pd.Series([i], index=[self._id])])

        # unmapped features will only contain it's own field
        unmapped_dict = {}
        for um in df_curated.index.difference(mapped_index):
            unmapped_dict[um] = {self._id: um}

        link.feature(
            dobject_id=dobject_id,
            features={**mapped_dict, **unmapped_dict},
            feature_entity=self._model.entity,
            species=self._model.species,
            featureset_name=self._featureset_name,
        )
