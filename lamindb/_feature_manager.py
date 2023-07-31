from collections import defaultdict
from typing import List, Optional, Union

import pandas as pd
from lamin_utils import logger
from lnschema_core.models import ORM, Dataset, Feature, FeatureSet, File

from ._queryset import QuerySet
from ._save import save


def validate_and_cast_feature(
    feature: Union[str, Feature], records: List[ORM]
) -> Feature:
    if isinstance(feature, str):
        feature_name = feature
        feature = Feature.filter(name=feature_name).one_or_none()
        if feature is None:
            orm_types = set([record.__class__ for record in records])
            feature = Feature(
                name=feature_name, type="category", registries=list(orm_types)
            )
            logger.warning(f"Created & saved: {feature}")
            feature.save()
    return feature


def create_features_df(
    file: File, feature_sets: List[FeatureSet], exclude: bool = True
):
    features = []
    for feature_set in feature_sets:
        if exclude:
            features_df = feature_set.features.exclude(registries__isnull=True).df()
        else:
            features_df = feature_set.features.df()
        slots = file.feature_sets.through.objects.filter(
            file=file, feature_set=feature_set
        ).list("slot")
        for slot in slots:
            features_df["slot"] = slot
            features.append(features_df)
    features_df = pd.concat(features)
    return features_df.sort_values(["slot", "registries"])


class FeatureManager:
    """Feature manager."""

    def __init__(self, host: Union[File, Dataset]):
        self._host = host
        self._compute_slots()

    def _compute_slots(self) -> None:
        slot_feature_sets = (
            self._feature_set_df_with_slots().reset_index().set_index("slot")["id"]
        )
        self._slots = {
            slot: self._host.feature_sets.get(id=i)
            for slot, i in slot_feature_sets.items()
        }

    def __repr__(self) -> str:
        self._compute_slots()
        if len(self._slots) > 0:
            msg = "slots:\n"
            for slot, feature_set in self._slots.items():
                msg += f"    {slot}: {feature_set}\n"
            return msg
        else:
            return "No linked features."

    def __getitem__(self, slot) -> QuerySet:
        id = (
            self._host.feature_sets.through.objects.filter(
                file_id=self._host.id, slot=slot
            )
            .one()
            .feature_set_id
        )
        accessor_by_orm = {
            field.related_model.__name__: field.name
            for field in self._host._meta.related_objects
        }
        accessor_by_orm["Feature"] = "features"
        feature_set = self._host.feature_sets.filter(id=id).one()
        orm_name = feature_set.ref_field.split(".")[1]
        return getattr(feature_set, accessor_by_orm[orm_name]).all()

    def _feature_set_df_with_slots(self) -> pd.DataFrame:
        """Return DataFrame."""
        df = self._host.feature_sets.df()
        df.insert(
            0,
            "slot",
            self._host.feature_sets.through.objects.filter(file_id=self._host.id)
            .df()
            .set_index("feature_set_id")
            .slot,
        )
        return df

    def add_labels(
        self, records: Union[ORM, List[ORM]], feature: Optional[Union[str, ORM]] = None
    ):
        """Add one or several labels and associate them with a feature."""
        if isinstance(records, str) or not isinstance(records, List):
            records = [records]
        if isinstance(records[0], str):  # type: ignore
            raise ValueError(
                "Please pass a record (an ORM object), not a string, e.g., via: label"
                f" = ln.Label(name='{records[0]}')"  # type: ignore
            )
        if self._host._state.adding:
            raise ValueError("Please save the file or dataset before adding a label!")
        feature = validate_and_cast_feature(feature, records)
        records_by_feature_orm = defaultdict(list)
        for record in records:
            if feature is None:
                error_msg = (
                    "Please pass feature: add_labels(labels, feature='myfeature')"
                )
                record_feature = feature
                if hasattr(record, "_feature"):
                    record_feature = record._feature
                if record_feature is None:
                    raise ValueError(error_msg)
                # TODO: refactor so that we don't call the following line
                # repeatedly for the same feature
                record_feature = validate_and_cast_feature(record_feature, [record])
            else:
                record_feature = feature
            records_by_feature_orm[
                (record_feature, record.__class__.__get_name_with_schema__())
            ].append(record)
        accessor_by_orm = {
            field.related_model.__get_name_with_schema__(): field.name
            for field in self._host._meta.related_objects
        }
        accessor_by_orm["core.Label"] = "labels"
        # ensure all labels are saved
        save(records)
        for (feature, orm_name), records in records_by_feature_orm.items():
            getattr(self._host, accessor_by_orm[orm_name]).add(
                *records, through_defaults={"feature_id": feature.id}
            )
        accessor_by_orm = {
            field.related_model.__name__: field.name
            for field in self._host._meta.related_objects
        }
        accessor_by_orm["Feature"] = "features"
        feature_set_links = self._host.feature_sets.through.objects.filter(
            file_id=self._host.id
        )
        feature_set_ids = [link.feature_set_id for link in feature_set_links.all()]
        # get all linked features of type Feature
        feature_sets = FeatureSet.filter(id__in=feature_set_ids).all()
        linked_features_by_slot = {
            feature_set_links.filter(feature_set_id=feature_set.id)
            .one()
            .slot: feature_set.features.all()
            for feature_set in feature_sets
            if "core.Feature" in feature_set.ref_field
        }
        for (feature, orm_name), records in records_by_feature_orm.items():
            feature = validate_and_cast_feature(feature, records)
            logger.info(f"Linking feature {feature.name} to {orm_name}")
            if feature.registries is None:
                feature.registries = orm_name
            elif orm_name not in feature.registries:
                feature.registries += f"|{orm_name}"
            feature.save()
            # check whether we have to update the feature set that manages labels
            # (Feature) to account for a new feature
            found_feature = False
            for _, linked_features in linked_features_by_slot.items():
                if feature in linked_features:
                    found_feature = True
            if not found_feature:
                if "ext" not in linked_features_by_slot:
                    logger.info("Creating feature_set for slot 'ext' (external)")
                    feature_set = FeatureSet([feature], modality="meta")
                    feature_set.save()
                    self.add_feature_set(feature_set, slot="ext")
                else:
                    feature_set = self._slots["ext"]
                    logger.info(
                        f"Linking feature {feature.name} to feature set {feature_set}"
                    )
                    feature_set.features.add(feature)
                    feature_set.n += 1
                    feature_set.save()
        self._compute_slots()

    def add_feature_set(self, feature_set: FeatureSet, slot: str):
        if self._host._state.adding:
            raise ValueError(
                "Please save the file or dataset before adding a feature set!"
            )
        feature_set.save()
        self._host.feature_sets.through(
            file=self._host, feature_set=feature_set, slot=slot
        ).save()
