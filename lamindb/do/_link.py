from typing import Optional  # noqa

import pandas as pd
from tabulate import tabulate  # type: ignore

from .._logger import colors, logger
from ..dev.db import insert
from ._query import query
from ._update import update


class LinkFeatureModel:
    def __init__(self, feature_model, featureset_name: str = None) -> None:
        self._feature_model = feature_model
        self._featureset_name = featureset_name

    @property
    def entity(self):
        """Correspond to the feature entity table."""
        return self._feature_model.entity

    @property
    def id_type(self):
        """Type of id used for curation."""
        return self._feature_model._id_field

    @property
    def species(self):
        """Species."""
        return self._feature_model.species

    @property
    def df(self):
        """Reference table."""
        return self._feature_model.df

    def curate(self, df: pd.DataFrame):
        if self.id_type in df.columns:
            return self._feature_model.curate(df=df, column=self.id_type)
        else:
            logger.warning(f"{self.id_type} column not found, using index as features.")
            return self._feature_model.curate(df=df, column=None)

    def ingest(self, dobject_id, df_curated):
        """Ingest features."""
        mapped_df = self.df.loc[df_curated.index[df_curated["__curated__"]]].copy()
        mapped_dict = {}
        for i, row in mapped_df.iterrows():
            mapped_dict[i] = pd.concat([row, pd.Series([i], index=[self.id_type])])

        link.feature(
            dobject_id=dobject_id,
            values=mapped_dict,
            feature_entity=self.entity,
            species=self.species,
            featureset_name=self._featureset_name,
        )


class link:
    """Link features and metadata to data."""

    @classmethod
    def feature_model(
        cls, df: pd.DataFrame, feature_model, featureset_name: str = None
    ):
        fm = LinkFeatureModel(feature_model, featureset_name=featureset_name)
        df_curated = fm.curate(df)
        n = df_curated["__curated__"].count()
        n_mapped = df_curated["__curated__"].sum()
        log = {
            "feature": fm.id_type,
            "n_mapped": n_mapped,
            "percent_mapped": round(n_mapped / n * 100, 1),
            "unmapped": df_curated.index[~df_curated["__curated__"]],
        }
        return (fm, df_curated), log

    @classmethod
    def feature(
        cls,
        dobject_id,
        values: dict,
        feature_entity: str,
        species: str,
        featureset_name: str = None,
    ):
        """Annotate genes."""
        species_id = insert.species(common_name=species)
        featureset_id = insert.features(
            features_dict=values,
            feature_entity=feature_entity,
            species=species,
            featureset_name=featureset_name,
        )

        # use the geneset_id and readout_id to create an entry in biometa
        biometa_id = insert.biometa(
            dobject_id=dobject_id,
            featureset_id=featureset_id,
        )

        logs = [[str(featureset_id), str(biometa_id), str(species_id)]]
        log_table = tabulate(
            logs,
            headers=[
                colors.green("geneset.id"),
                colors.purple("biometa.id"),
                colors.blue("species.id"),
            ],
            tablefmt="pretty",
        )
        logger.success(
            f"Annotated data {dobject_id} with the following features:\n{log_table}",
        )

    @classmethod
    def readout(cls, dobject_id, efo_id: str):
        """Link readout."""
        readout_id = insert.readout(efo_id=efo_id)

        # query biometa associated with a dobject
        query_dobject_biometa = getattr(query, "dobject_biometa")
        dobject_biometa = query_dobject_biometa(dobject_id=dobject_id)
        if len(dobject_biometa) > 0:
            biometa_ids = [i.biometa_id for i in dobject_biometa]
        else:
            biometa_ids = [insert.biometa(dobject_id=dobject_id)]
            logger.warning(
                f"No biometa found for dobject {dobject_id}, created biometa"
                f" {biometa_ids[0]}"
            )

        # fill in biometa entries with readout_id
        for biometa_id in biometa_ids:
            update_biometa = getattr(update, "biometa")
            update_biometa(biometa_id, readout_id=readout_id)

        logger.success(
            f"{colors.blue(f'readout_id {readout_id}')} has been added to"
            f" {colors.purple(f'biometa entries {biometa_ids}')} associated with"
            f" {colors.green(f'dobject {dobject_id}')}."
        )
