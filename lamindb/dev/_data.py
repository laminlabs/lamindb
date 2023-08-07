from collections import defaultdict
from typing import Dict, List, Optional, Union

from lamin_utils import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core.models import Data, Feature, FeatureSet, Label, Registry

from .._query_set import QuerySet
from .._registry import get_default_str_field
from .._save import save
from ._feature_manager import FeatureManager
from .exc import ValidationError


def validate_and_cast_feature(
    feature: Union[str, Feature], records: List[Registry]
) -> Feature:
    if isinstance(feature, str):
        feature_name = feature
        feature = Feature.filter(name=feature_name).one_or_none()
        if feature is None:
            registries = set(
                [record.__class__.__get_name_with_schema__() for record in records]
            )
            registries_str = "|".join(registries)
            msg = (
                f"ln.Feature(name='{feature_name}', type='category',"
                f" registries='{registries_str}').save()"
            )
            raise ValidationError(f"Feature not validated. If it looks correct: {msg}")
    return feature


@doc_args(Data.get_labels.__doc__)
def get_labels(
    self,
    feature: Optional[Union[str, Registry]] = None,
    mute: bool = False,
    flat_names: bool = False,
) -> Union[QuerySet, Dict[str, QuerySet], List]:
    """{}"""
    if isinstance(feature, str):
        feature_name = feature
        feature = Feature.filter(name=feature_name).one_or_none()
        if feature is None:
            raise ValueError("feature doesn't exist")
    if feature.registries is None:
        raise ValueError("feature does not have linked labels")
    registries_to_check = feature.registries.split("|")
    if len(registries_to_check) > 1 and not mute:
        logger.warning("labels come from multiple registries!")
    qs_by_registry = {}
    for registry in registries_to_check:
        # currently need to distinguish between Label and non-Label, because
        # we only have the feature information for Label
        if registry == "core.Label":
            links_to_labels = getattr(
                self, self.features._accessor_by_orm[registry]
            ).through.objects.filter(file_id=self.id, feature_id=feature.id)
            label_ids = [link.label_id for link in links_to_labels]
            qs_by_registry[registry] = Label.objects.filter(id__in=label_ids)
        else:
            qs_by_registry[registry] = getattr(
                self, self.features._accessor_by_orm[registry]
            ).all()
    if flat_names:
        # returns a flat list of names
        from .._registry import get_default_str_field

        values = []
        for v in qs_by_registry.values():
            values += v.list(get_default_str_field(v))
        return values
    if len(registries_to_check) == 1:
        return qs_by_registry[registry]
    else:
        return qs_by_registry


@doc_args(Data.add_labels.__doc__)
def add_labels(
    self,
    records: Union[Registry, List[Registry], QuerySet],
    feature: Optional[Union[str, Registry]] = None,
) -> None:
    """{}"""
    if isinstance(records, (QuerySet, QuerySet.__base__)):  # need to have both
        records = records.list()
    if isinstance(records, str) or not isinstance(records, List):
        records = [records]
    if isinstance(records[0], str):  # type: ignore
        raise ValueError(
            "Please pass a record (a `Registry` object), not a string, e.g., via:"
            " label"
            f" = ln.Label(name='{records[0]}')"  # type: ignore
        )
    if self._state.adding:
        raise ValueError("Please save the file or dataset before adding a label!")
    for record in records:
        if record._state.adding:
            raise ValidationError(
                f"{record} not validated. If it looks correct: record.save()"
            )
    feature = validate_and_cast_feature(feature, records)
    orig_feature = feature
    records_by_feature_orm = defaultdict(list)
    for record in records:
        if feature is None:
            error_msg = "Please pass feature: add_labels(labels, feature='myfeature')"
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
        getattr(self, self.features._accessor_by_orm[orm_name]).add(
            *records, through_defaults={"feature_id": feature.id}
        )
    feature_set_links = self.feature_sets.through.objects.filter(file_id=self.id)
    feature_set_ids = [link.feature_set_id for link in feature_set_links.all()]
    # get all linked features of type Feature
    feature_sets = FeatureSet.filter(id__in=feature_set_ids).all()
    linked_features_by_slot = {
        feature_set_links.filter(feature_set_id=feature_set.id)
        .one()
        .slot: feature_set.features.all()
        for feature_set in feature_sets
        if "core.Feature" == feature_set.registry
    }
    for (feature, orm_name), records in records_by_feature_orm.items():
        feature = validate_and_cast_feature(feature, records)
        msg = ""
        if orig_feature is None:
            records_display = ", ".join(
                [
                    f"'{getattr(record, get_default_str_field(record))}'"
                    for record in records
                ]
            )
            msg += f"linked labels {records_display} to feature '{feature.name}'"
        if feature.registries is None or orm_name not in feature.registries:
            if len(msg) > 0:
                msg += ", "
            msg += f"linked feature '{feature.name}' to registry '{orm_name}'"
            if feature.registries is None:
                feature.registries = orm_name
            elif orm_name not in feature.registries:
                feature.registries += f"|{orm_name}"
            feature.save()
        if len(msg) > 0:
            logger.save(msg)
        # check whether we have to update the feature set that manages labels
        # (Feature) to account for a new feature
        found_feature = False
        for _, linked_features in linked_features_by_slot.items():
            if feature in linked_features:
                found_feature = True
        if not found_feature:
            if "external" not in linked_features_by_slot:
                feature_set = FeatureSet([feature], modality="meta")
                feature_set.save()
                self.features.add_feature_set(feature_set, slot="external")
                logger.save("created feature set for slot 'external'")
            else:
                feature_set = self.features._feature_set_by_slot["external"]
                feature_set.features.add(feature)
                feature_set.n += 1
                feature_set.save()
                logger.save(
                    f"linked feature {feature.name} to feature set {feature_set}"
                )


@property  # type: ignore
@doc_args(Data.features.__doc__)
def features(self) -> "FeatureManager":
    """{}"""
    from lamindb.dev._feature_manager import FeatureManager

    return FeatureManager(self)


setattr(Data, "features", features)
setattr(Data, "add_labels", add_labels)
setattr(Data, "get_labels", get_labels)
