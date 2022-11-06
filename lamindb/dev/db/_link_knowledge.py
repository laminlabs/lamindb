from typing import List, Union

import lnschema_core as core
import pandas as pd
from lamin_logger import colors, logger
from tabulate import tabulate  # type: ignore

from lamindb._link import link
from lamindb.dev.db._add import add
from lamindb.dev.db._select import select
from lamindb.knowledge import CellMarker, Gene, Protein
from lamindb.schema import bionty, wetlab
from lamindb.schema._table import table_meta


class LinkFeatureToKnowledgeTable:
    """Link features to a knowledge table.

    Args:
        df: a DataFrame
        knowledge_table: knowledge table
        featureset_name: name of the featureset
            If provided, will create a featureset if the same name does not exist
    """

    def __init__(
        self,
        df: pd.DataFrame,
        knowledge_table: Union[Gene, Protein, CellMarker],
        featureset_name: str = None,
    ) -> None:
        self._knowledge = knowledge_table
        self._name = self._knowledge.__class__.__name__
        self._id = knowledge_table._id_field
        self._featureset_name = featureset_name

        # Add and curate features against a knowledge table
        column = None
        if self._id in df.columns:
            column = self._id
        else:
            logger.warning(f"{self._id} column not found, using index as features.")
        df_curated = self._knowledge.curate(df=df, column=column)

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

        # create a featureset record if featureset_name is provided
        self._featureset = None
        if self._featureset_name is not None:
            self._featureset = select(
                bionty.Featureset,
                name=self._featureset_name,
                feature_entity=self._knowledge.entity,
            ).one_or_none()

    @property
    def log(self) -> dict:
        """Logging."""
        return self._log

    def commit(self, dobject: core.DObject) -> None:
        """Commit features."""
        # insert species entry if not exists
        species = select(
            bionty.Species, common_name=self._knowledge.species
        ).one_or_none()
        if species is None:
            species = add(bionty.Species(common_name=self._knowledge.species))

        model = table_meta.get_model(self._name)

        # all existing feature rows of a species in the db
        db_rows = (
            select(model)
            .where(getattr(model, self._id).in_(self._df_curated.index))
            .where(getattr(model, "species_id") == species.id)
            .df()
        )

        # ids of the existing features
        exist_features = self._df_curated.index.intersection(db_rows[self._id])
        feature_ids = set(db_rows[db_rows[self._id].isin(exist_features)].index)

        # new features to be inserted
        new_ids = self._df_curated.index.difference(exist_features)
        records = []
        if len(new_ids) > 0:
            # mapped new_ids
            mapped = self._knowledge.df.loc[
                self._knowledge.df.index.intersection(new_ids)
            ].copy()
            mapped.index.name = self._id
            if mapped.shape[0] > 0:
                for kwargs in mapped.reset_index().to_dict(orient="records"):
                    record = getattr(bionty, f"{self._name}")(**kwargs)
                    records.append(record)
                    feature_ids.add(record.id)
            # unmapped new_ids
            unmapped = set(new_ids).difference(mapped.index)
            if len(unmapped) > 0:
                for i in unmapped:
                    record = getattr(bionty, f"{self._name}")(
                        **{self._id: i, "species_id": species.id}
                    )
                    records.append(record)
                    feature_ids.add(record.id)
        else:
            # search if a featureset containing the same genes already exists
            featuresets = select(getattr(bionty, f"Featureset{self._name}")).df()
            for i, fset in (
                featuresets.groupby("featureset_id", group_keys=True)[
                    f"{self._knowledge.entity}_id"
                ]
                .apply(set)
                .items()
            ):
                if fset == feature_ids:
                    self._featureset = select(bionty.Featureset, id=i).one()
                    break

        # 1) insert the featureset
        # if a featureset already exist, no features or link entries will be created
        if self._featureset is None:
            self._featureset = add(
                bionty.Featureset(
                    name=self._featureset_name, feature_entity=self._knowledge.entity
                )
            )
            # 2) insert the feature records
            if len(records) > 0:
                add(records)

            # 3) insert rows into the link table
            link_records = [
                getattr(bionty, f"Featureset{self._name}")(
                    **{
                        "featureset_id": self._featureset.id,
                        f"{self._knowledge.entity}_id": i,
                    }
                )
                for i in feature_ids
            ]
            add(link_records)

        # 4) use the featureset_id to create an entry in biometa
        dobject_biometas = select(wetlab.DObjectBiometa, dobject_id=dobject.id).all()
        if len(dobject_biometas) == 0:
            biometa = add(wetlab.Biometa(featureset_id=self._featureset.id))
            link(dobject, biometa)
        else:
            raise NotImplementedError

        logs: List = [[str(self._featureset.id), str(biometa.id), species.id]]
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
            f"Annotated data {dobject.id} with the following features:\n{log_table}",
        )
