from typing import Dict, List, Union

import pandas as pd
from lamin_utils import logger
from lnschema_core.models import Dataset, Feature, FeatureSet, File

from .._query_set import QuerySet


def get_host_id_field(host: Union[File, Dataset]) -> str:
    if isinstance(host, File):
        host_id_field = "file_id"
    else:
        host_id_field = "dataset_id"
    return host_id_field


def create_features_df(
    host: Union[File, Dataset], feature_sets: List[FeatureSet], exclude: bool = True
):
    features = []
    for feature_set in feature_sets:
        if exclude:
            features_df = feature_set.features.exclude(registries__isnull=True).df()
        else:
            features_df = feature_set.features.df()
        host_id_field = get_host_id_field(host)
        kwargs = {host_id_field: host.id, "feature_set": feature_set}
        slots = host.feature_sets.through.objects.filter(**kwargs).list("slot")
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
    # if the host is not yet saved
    if host._state.adding:
        return host._feature_sets
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id}
    # otherwise, we need a query
    feature_set_links = host.feature_sets.through.objects.filter(**kwargs)
    return {
        feature_set_link.slot: FeatureSet.objects.get(
            id=feature_set_link.feature_set_id
        )
        for feature_set_link in feature_set_links
    }


def get_label_links(
    host: Union[File, Dataset], registry: str, feature: Feature
) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id, "feature_id": feature.id}
    link_records = getattr(
        host, host.features._accessor_by_orm[registry]
    ).through.objects.filter(**kwargs)
    return link_records


def get_feature_set_links(host: Union[File, Dataset]) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id}
    feature_set_links = host.feature_sets.through.objects.filter(**kwargs)
    return feature_set_links


class FeatureManager:
    """Feature manager (:attr:`~lamindb.dev.Data.features`).

    See :class:`~lamindb.dev.Data` for more information.
    """

    def __init__(self, host: Union[File, Dataset]):
        self._host = host
        self._feature_set_by_slot = get_feature_set_by_slot(host)
        self._accessor_by_orm = get_accessor_by_orm(host)

    def __repr__(self) -> str:
        if len(self._feature_set_by_slot) > 0:
            msg = ""
            no_modality = []
            for slot, feature_set in self._feature_set_by_slot.items():
                msg += f"'{slot}': {feature_set}\n"
                if feature_set.modality is None:
                    no_modality.append(feature_set.id)
            if len(no_modality) > 0:
                ids = ", ".join(f"'{key}'" for key in no_modality)
                s = "" if len(no_modality) == 1 else "s"
                logger.debug(f"consider assigning modality to feature set{s}: {ids}")
            return msg
        else:
            return "no linked features"

    def __getitem__(self, slot) -> QuerySet:
        if slot not in self._feature_set_by_slot:
            raise ValueError(
                f"No linked feature set for slot: {slot}\nDid you get validation"
                " warnings? Only features that match registered features get validated"
                " and linked."
            )
        feature_set = self._feature_set_by_slot[slot]
        orm_name = feature_set.registry
        return getattr(feature_set, self._accessor_by_orm[orm_name]).all()

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
        host_id_field = get_host_id_field(self._host)
        kwargs = {
            host_id_field: self._host.id,
            "feature_set": feature_set,
            "slot": slot,
        }
        link_record = self._host.feature_sets.through.objects.filter(
            **kwargs
        ).one_or_none()
        if link_record is None:
            self._host.feature_sets.through(**kwargs).save()
            self._feature_set_by_slot[slot] = feature_set

    def get_feature_set(self, slot: str) -> FeatureSet:
        """Get feature set by slot.

        Args:
            slot: `str` The access slot.
        """
        return self._feature_set_by_slot.get(slot)
