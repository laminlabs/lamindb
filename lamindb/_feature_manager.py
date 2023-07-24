from collections import defaultdict
from typing import List, Union

from lnschema_core.models import ORM, Dataset, File

from ._save import save


class FeatureManager:
    """Feature manager."""

    def __init__(self, host: Union[File, Dataset]):
        self._host = host

    def add_labels(self, records: List[ORM]):
        """Add new labels and associate them with a feature."""
        records_by_orm = defaultdict(list)
        records_by_feature_orm = defaultdict(list)
        for record in records:
            records_by_orm[record.__class__.__name__].append(record)
            feature = record._feature if hasattr(record, "_feature") else record.feature
            records_by_feature_orm[(feature, record.__class__.__name__)].append(record)
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
            feature.labels_orm = orm_name
            feature.labels_schema = schema_and_accessor_by_orm[orm_name][0]
            feature.save()
