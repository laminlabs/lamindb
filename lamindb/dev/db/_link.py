import pandas as pd
from lamin_logger import colors, logger
from tabulate import tabulate  # type: ignore
from typing_extensions import Literal

from lamindb.dev._features import add_features_and_featureset
from lamindb.dev.db._add import add
from lamindb.dev.db._select import select
from lamindb.schema import bionty, wetlab


class LinkFeatureModel:
    """Link a feature model to dobject during ingestion.

    Args:
        knowledge_table: a feature model instance
        featureset_name: name of the featureset

    For an overview of feature models, see: `bionty.lookup.knowledge_table <https://lamin.ai/docs/bionty/bionty.lookup#bionty.lookup.knowledge_table)>`__.
    """  # noqa

    def __init__(self, knowledge_table, featureset_name: str = None) -> None:
        self._knowledge_table = knowledge_table
        self._featureset_name = featureset_name

    @property
    def entity(self):
        """Correspond to the feature entity table."""
        return self._knowledge_table.entity

    @property
    def id_type(self):
        """Type of id used for curation."""
        return self._knowledge_table._id_field

    @property
    def species(self):
        """Species."""
        return self._knowledge_table.species

    @property
    def df(self):
        """Reference table."""
        return self._knowledge_table.df

    def curate(self, df: pd.DataFrame):
        column = None
        if self.id_type in df.columns:
            column = self.id_type
        else:
            logger.warning(f"{self.id_type} column not found, using index as features.")
        df_curated = self._knowledge_table.curate(df=df, column=column)

        # logging of curation
        n = df_curated["__curated__"].count()
        n_mapped = df_curated["__curated__"].sum()
        log = {
            "feature": self.id_type,
            "n_mapped": n_mapped,
            "percent_mapped": round(n_mapped / n * 100, 1),
            "unmapped": df_curated.index[~df_curated["__curated__"]],
        }

        return df_curated, log

    def ingest(self, dobject_id: str, df_curated: pd.DataFrame):
        """Ingest features."""
        # mapped features will also contain fields in the reference table
        mapped_index = df_curated.index[df_curated["__curated__"]]
        mapped_df = self.df.loc[mapped_index].copy()
        mapped_dict = {}
        for i, row in mapped_df.iterrows():
            mapped_dict[i] = pd.concat([row, pd.Series([i], index=[self.id_type])])

        # unmapped features will only contain it's own field
        unmapped_dict = {}
        for um in df_curated.index.difference(mapped_index):
            unmapped_dict[um] = {self.id_type: um}

        link.feature(
            dobject_id=dobject_id,
            features={**mapped_dict, **unmapped_dict},
            feature_entity=self.entity,
            species=self.species,
            featureset_name=self._featureset_name,
        )


class link:
    """Link metadata.

    Guide: :doc:`/db/guide/link`.
    """

    @classmethod
    def knowledge_table(
        cls, df: pd.DataFrame, knowledge_table, featureset_name: str = None
    ):
        """Curate a DataFrame with a feature model."""
        fm = LinkFeatureModel(knowledge_table, featureset_name=featureset_name)
        df_curated, log = fm.curate(df)

        return {"model": fm, "df_curated": df_curated, "log": log}

    @classmethod
    def feature(
        cls,
        dobject_id: str,
        features: dict,  # what is a features dict? can we have something more typed?
        feature_entity: Literal["gene", "protein", "cell_marker"],
        species: bionty.species,
        featureset_name: str = None,
    ):
        """Annotate dobject with features.

        Add all features and a featureset.

        Link featureset to biometa.

        Link biometa to dobject.
        """
        featureset = add_features_and_featureset(
            features=features,
            feature_entity=feature_entity,
            species=species.common_name,
            name=featureset_name,
        )

        # use the featureset_id to create an entry in biometa
        # TODO: need to make this easier
        dobject_biometas = select(wetlab.dobject_biometa, dobject_id=dobject_id).all()
        if len(dobject_biometas) == 0:
            biometa = add(bionty.biometa(featureset_id=featureset.id))
            cls.biometa(dobject_id=dobject_id, biometa_id=biometa.id)
        else:
            raise NotImplementedError

        logs = [[str(featureset.id), str(biometa.id), str(species.id)]]  # type: ignore
        log_table = tabulate(
            logs,
            headers=[
                colors.green("featureset.id"),
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
        readout = add(wetlab.readout(efo_id=efo_id))

        # select biometa associated with a dobject
        dobject_biometa = select(wetlab.dobject_biometa, dobject_id=dobject_id).all()
        if len(dobject_biometa) > 0:
            biometa_ids = [i.biometa_id for i in dobject_biometa]
        else:
            # TODO: fix here
            biometa_ids = [add(wetlab.biometa(readout_id=readout.id)).id]
            logger.warning(
                f"No biometa found for dobject {dobject_id}, created biometa"
                f" {biometa_ids[0]}"
            )

        # fill in biometa entries with readout_id
        for biometa_id in biometa_ids:
            add(wetlab.biometa(biometa_id, readout_id=readout.id))

        logger.success(
            f"Added {colors.blue(f'readout_id {readout.id}')} to"
            f" {colors.purple(f'biometa {biometa_ids}')} linked to"
            f" {colors.green(f'dobject {dobject_id}')}."
        )

    @classmethod
    def biometa(cls, dobject_id: str, biometa_id: str):
        """Link a dobject to a biometa."""
        dobject_biometas = select(
            wetlab.dobject_biometa, dobject_id=dobject_id, biometa_id=biometa_id
        ).all()
        if len(dobject_biometas) > 0:
            raise AssertionError(
                f"dobject {dobject_id} is already linked to biometa {biometa_id}!"
            )
        else:
            add(wetlab.dobject_biometa(dobject_id=dobject_id, biometa_id=biometa_id))
