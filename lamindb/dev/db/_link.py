import pandas as pd
from lamin_logger import colors, logger
from tabulate import tabulate  # type: ignore

from ._insert import insert
from ._select import select
from ._update import update


class LinkFeatureModel:
    """Link a feature model to dobject during ingestion.

    Args:
        feature_model: a feature model instance
        featureset_name: name of the featureset

    For an overview of feature models, see: `bionty.lookup.feature_model <https://lamin.ai/docs/bionty/bionty.lookup#bionty.lookup.feature_model)>`__.
    """  # noqa

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
        column = None
        if self.id_type in df.columns:
            column = self.id_type
        else:
            logger.warning(f"{self.id_type} column not found, using index as features.")
        df_curated = self._feature_model.curate(df=df, column=column)

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
            values={**mapped_dict, **unmapped_dict},
            feature_entity=self.entity,
            species=self.species,
            featureset_name=self._featureset_name,
        )


class link:
    """Link metadata.

    Guide: :doc:`/db/guide/link`.
    """

    @classmethod
    def feature_model(
        cls, df: pd.DataFrame, feature_model, featureset_name: str = None
    ):
        """Curate a DataFrame with a feature model."""
        fm = LinkFeatureModel(feature_model, featureset_name=featureset_name)
        df_curated, log = fm.curate(df)

        return {"model": fm, "df_curated": df_curated, "log": log}

    @classmethod
    def feature(
        cls,
        dobject_id: str,
        values: dict,
        feature_entity: str,
        species: str,
        featureset_name: str = None,
    ):
        """Annotate genes."""
        species = insert.species(common_name=species)  # type: ignore

        featureset = insert.featureset_from_features(  # type: ignore
            features_dict=values,
            feature_entity=feature_entity,
            species=species.common_name,  # type: ignore
            featureset_name=featureset_name,
        )

        # use the featureset_id to create an entry in biometa
        # TODO: need to make this easier
        dobject_biometas = select.dobject_biometa(dobject_id=dobject_id).all()  # type: ignore  # noqa
        if len(dobject_biometas) == 0:
            # insert a biometa entry and link to dobject
            # TODO: force insert here
            biometa = getattr(insert, "biometa")(
                featureset_id=featureset.id, force=True
            )
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
        readout = insert.readout(efo_id=efo_id)  # type: ignore

        # select biometa associated with a dobject
        dobject_biometa = select.dobject_biometa(dobject_id=dobject_id).all()  # type: ignore  # noqa
        if len(dobject_biometa) > 0:
            biometa_ids = [i.biometa_id for i in dobject_biometa]
        else:
            # TODO: fix here
            biometa_ids = [insert.biometa(readout_id=readout.id).id]  # type: ignore
            logger.warning(
                f"No biometa found for dobject {dobject_id}, created biometa"
                f" {biometa_ids[0]}"
            )

        # fill in biometa entries with readout_id
        for biometa_id in biometa_ids:
            update_biometa = getattr(update, "biometa")
            update_biometa(biometa_id, readout_id=readout.id)

        logger.success(
            f"Added {colors.blue(f'readout_id {readout.id}')} to"
            f" {colors.purple(f'biometa {biometa_ids}')} linked to"
            f" {colors.green(f'dobject {dobject_id}')}."
        )

    @classmethod
    def biometa(cls, dobject_id: str, biometa_id: int):
        """Link a dobject to a biometa."""
        dobject_biometas = getattr(select, "dobject_biometa")(
            dobject_id=dobject_id, biometa_id=biometa_id
        ).all()
        if len(dobject_biometas) > 0:
            raise AssertionError(
                f"dobject {dobject_id} is already linked to biometa {biometa_id}!"
            )
        else:
            _ = getattr(insert, "dobject_biometa")(
                dobject_id=dobject_id, biometa_id=biometa_id
            )
