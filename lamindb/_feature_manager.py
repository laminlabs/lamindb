from collections import defaultdict
from typing import Dict, List, Optional, Union

import pandas as pd
from lamin_utils import logger
from lnschema_core.models import ORM, Dataset, Feature, FeatureSet, File, Label

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


def get_accessor_by_orm(host: Union[File, Dataset]) -> Dict:
    dictionary = {
        field.related_model.__get_name_with_schema__(): field.name
        for field in host._meta.related_objects
    }
    dictionary["core.Feature"] = "features"
    dictionary["core.Label"] = "labels"
    return dictionary


def get_feature_set_by_slot(host) -> Dict:
    feature_set_links = host.feature_sets.through.objects.filter(file_id=host.id)
    return {
        feature_set_link.slot: FeatureSet.objects.get(
            id=feature_set_link.feature_set_id
        )
        for feature_set_link in feature_set_links
    }


class FeatureManager:
    """Feature manager."""

    def __init__(self, host: Union[File, Dataset]):
        self._host = host
        self._feature_set_by_slot = get_feature_set_by_slot(host)
        self._accessor_by_orm = get_accessor_by_orm(host)

    def __repr__(self) -> str:
        if len(self._feature_set_by_slot) > 0:
            msg = "slots:\n"
            for slot, feature_set in self._feature_set_by_slot.items():
                msg += f"    {slot}: {feature_set}\n"
            return msg
        else:
            return "No linked features."

    def __getitem__(self, slot) -> QuerySet:
        feature_set = self._feature_set_by_slot[slot]
        orm_name = ".".join(feature_set.ref_field.split(".")[:2])
        return getattr(feature_set, self._accessor_by_orm[orm_name]).all()

    def get_labels(
        self,
        feature: Optional[Union[str, ORM]] = None,
        mute: bool = False,
        flat_names: bool = False,
    ) -> Union[QuerySet, Dict[str, QuerySet], List]:
        """Get labels given a feature."""
        if isinstance(feature, str):
            feature_name = feature
            feature = Feature.filter(name=feature_name).one_or_none()
            if feature is None:
                raise ValueError("Feature doesn't exist")
        if feature.registries is None:
            raise ValueError("Feature does not have linked labels")
        registries_to_check = feature.registries.split("|")
        if len(registries_to_check) > 1 and not mute:
            logger.warning("Labels come from multiple registries!")
        qs_by_registry = {}
        for registry in registries_to_check:
            # currently need to distinguish between Label and non-Label, because
            # we only have the feature information for Label
            if registry == "core.Label":
                links_to_labels = getattr(
                    self._host, self._accessor_by_orm[registry]
                ).through.objects.filter(file_id=self._host.id, feature_id=feature.id)
                label_ids = [link.label_id for link in links_to_labels]
                qs_by_registry[registry] = Label.objects.filter(id__in=label_ids)
            else:
                qs_by_registry[registry] = getattr(
                    self._host, self._accessor_by_orm[registry]
                ).all()
        if flat_names:
            # returns a flat list of names
            from ._orm import get_default_str_field

            values = []
            for v in qs_by_registry.values():
                values += v.list(get_default_str_field(v))
            return values
        if len(registries_to_check) == 1:
            return qs_by_registry[registry]
        else:
            return qs_by_registry

    def add_labels(
        self,
        records: Union[ORM, List[ORM], QuerySet],
        feature: Optional[Union[str, ORM]] = None,
    ) -> None:
        """Add one or several labels and associate them with a feature."""
        if isinstance(records, (QuerySet, QuerySet.__base__)):  # need to have both
            records = records.list()
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
        # ensure all labels are saved
        save(records)
        for (feature, orm_name), records in records_by_feature_orm.items():
            getattr(self._host, self._accessor_by_orm[orm_name]).add(
                *records, through_defaults={"feature_id": feature.id}
            )
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
                    feature_set = self._feature_set_by_slot["ext"]
                    logger.info(
                        f"Linking feature {feature.name} to feature set {feature_set}"
                    )
                    feature_set.features.add(feature)
                    feature_set.n += 1
                    feature_set.save()

    def add_feature_set(self, feature_set: FeatureSet, slot: str):
        """Add new feature set to a slot.

        Args:
            feature_set: `FeatureSet` A feature set object.
            slot: `str` The access slot.
        """
        if self._host._state.adding:
            raise ValueError(
                "Please save the file or dataset before adding a feature set!"
            )
        feature_set.save()
        self._host.feature_sets.through(
            file=self._host, feature_set=feature_set, slot=slot
        ).save()
        self._feature_set_by_slot[slot] = feature_set
