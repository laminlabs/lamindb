from collections import defaultdict
from typing import List, Optional, Union

import pandas as pd
from lamin_utils import logger
from lnschema_core.models import ORM, Dataset, Feature, File

from ._queryset import QuerySet
from ._save import save


class FeatureManager:
    """Feature manager."""

    def __init__(self, host: Union[File, Dataset]):
        self._host = host

    def add_labels(self, records: List[ORM], feature: Optional[Union[str, ORM]] = None):
        """Add new labels and associate them with a feature."""
        if isinstance(feature, str):
            feature = Feature.select(name=feature).one()
        records_by_orm = defaultdict(list)
        records_by_feature_orm = defaultdict(list)
        for record in records:
            records_by_orm[record.__class__.__name__].append(record)
            if feature is None:
                try:
                    record_feature = (
                        record._feature
                        if hasattr(record, "_feature")
                        else record.feature
                    )
                except ValueError:
                    raise ValueError("Pass feature argument")
            else:
                record_feature = feature
            records_by_feature_orm[(record_feature, record.__class__.__name__)].append(
                record
            )
        schema_and_accessor_by_orm = {
            field.related_model.__name__: (
                field.related_model.__get_schema_name__(),
                field.name,
            )
            for field in self._host._meta.related_objects
        }
        schema_and_accessor_by_orm["Label"] = ("core", "labels")
        for orm_name, records in records_by_orm.items():
            save(records)
            getattr(self._host, schema_and_accessor_by_orm[orm_name][1]).set(records)
        for (feature, orm_name), records in records_by_feature_orm.items():
            logger.info(f"Linking feature {feature.name} to {orm_name}")
            feature.labels_orm = orm_name
            feature.labels_schema = schema_and_accessor_by_orm[orm_name][0]
            feature.save()

    def all(self, slot: str = "columns") -> QuerySet:
        """All features for a slot."""
        id = (
            self._host.feature_sets.through.objects.filter(
                file_id=self._host.id, slot=slot
            )
            .one()
            .featureset_id
        )
        accessor_by_orm = {
            field.related_model.__name__: field.name
            for field in self._host._meta.related_objects
        }
        accessor_by_orm["Feature"] = "features"
        feature_set = self._host.feature_sets.filter(id=id).one()
        return getattr(feature_set, accessor_by_orm[feature_set.ref_orm]).all()

    def df(self) -> pd.DataFrame:
        """Return DataFrame."""
        df = self._host.feature_sets.df()
        df.insert(
            0,
            "slot",
            self._host.feature_sets.through.objects.filter(file_id=self._host.id)
            .df()
            .set_index("featureset_id")
            .slot,
        )
        return df
