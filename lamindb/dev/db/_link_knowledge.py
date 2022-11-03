from typing import List, Union

import pandas as pd
from lamin_logger import colors, logger
from tabulate import tabulate  # type: ignore

from lamindb._link import link
from lamindb.dev.db._add import add
from lamindb.dev.db._select import select
from lamindb.knowledge import CellMarker, Gene, Protein, Species
from lamindb.schema import bionty, core, wetlab
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

        # create a featureset record if featureset_name is provided
        self._featureset = None
        if self._featureset_name is not None:
            self._featureset = select(
                bionty.featureset,
                name=self._featureset_name,
                feature_entity=self._model.entity,
            ).one_or_none()

    @property
    def log(self) -> dict:
        """Logging."""
        return self._log

    def commit(self, dobject: core.DObject) -> None:
        """Commit features."""
        # insert species entry if not exists
        # TODO: insert species
        species = select(bionty.species, common_name=self._model.species).one_or_none()
        if species is None:
            kwargs = Species().df.loc[self._model.species].to_dict()
            kwargs["common_name"] = self._model.species
            species = add(bionty.species(**kwargs))

        model = table_meta.get_model(self._model.entity)

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
            mapped = self._model.df.loc[
                self._model.df.index.intersection(new_ids)
            ].copy()
            mapped.index.name = self._id
            if mapped.shape[0] > 0:
                for kwargs in mapped.reset_index().to_dict(orient="records"):
                    record = getattr(bionty, f"{self._model.entity}")(**kwargs)
                    records.append(record)
                    feature_ids.add(record.id)
            # unmapped new_ids
            unmapped = set(new_ids).difference(mapped.index)
            if len(unmapped) > 0:
                for i in unmapped:
                    record = getattr(bionty, f"{self._model.entity}")(
                        **{self._id: i, "species_id": species.id}
                    )
                    records.append(record)
                    feature_ids.add(record.id)
        else:
            # search if a featureset containing the same genes already exists
            featuresets = select(
                getattr(bionty, f"featureset_{self._model.entity}")
            ).df()
            for i, fset in (
                featuresets.groupby("featureset_id", group_keys=True)[
                    f"{self._model.entity}_id"
                ]
                .apply(set)
                .items()
            ):
                if fset == feature_ids:
                    self._featureset = select(bionty.featureset, id=i).one()
                    break

        # 1) insert the featureset
        # if a featureset already exist, no features or link entries will be created
        if self._featureset is None:
            self._featureset = add(
                bionty.featureset(
                    name=self._featureset_name, feature_entity=self._model.entity
                )
            )
            # 2) insert the feature records
            if len(records) > 0:
                add(records)

            # 3) insert rows into the link table
            link_records = [
                getattr(bionty, f"featureset_{self._model.entity}")(
                    **{
                        "featureset_id": self._featureset.id,
                        f"{self._model.entity}_id": i,
                    }
                )
                for i in feature_ids
            ]
            add(link_records)

        # 4) use the featureset_id to create an entry in biometa
        dobject_biometas = select(wetlab.dobject_biometa, dobject_id=dobject.id).all()
        if len(dobject_biometas) == 0:
            biometa = add(wetlab.biometa(featureset_id=self._featureset.id))
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
