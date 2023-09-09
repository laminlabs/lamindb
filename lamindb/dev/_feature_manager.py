from typing import Dict, Union

from lamin_utils import colors
from lnschema_core.models import Data, Dataset, Feature, FeatureSet, File

from .._query_set import QuerySet
from .._registry import get_default_str_field


def dict_related_model_to_related_name(orm):
    d: Dict = {
        i.related_model.__get_name_with_schema__(): i.related_name
        for i in orm._meta.related_objects
        if i.related_name is not None
    }
    d.update(
        {
            i.related_model.__get_name_with_schema__(): i.name
            for i in orm._meta.many_to_many
            if i.name is not None
        }
    )

    return d


def dict_schema_name_to_model_name(orm):
    d: Dict = {
        i.related_model.__get_name_with_schema__(): i.related_model
        for i in orm._meta.related_objects
        if i.related_name is not None
    }
    d.update(
        {
            i.related_model.__get_name_with_schema__(): i.related_model
            for i in orm._meta.many_to_many
            if i.name is not None
        }
    )

    return d


def get_host_id_field(host: Union[File, Dataset]) -> str:
    if isinstance(host, File):
        host_id_field = "file_id"
    else:
        host_id_field = "dataset_id"
    return host_id_field


def get_accessor_by_orm(host: Union[File, Dataset]) -> Dict:
    dictionary = {
        field.related_model.__get_name_with_schema__(): field.name
        for field in host._meta.related_objects
    }
    dictionary["core.Feature"] = "features"
    dictionary["core.ULabel"] = "ulabels"
    return dictionary


def get_feature_set_by_slot(host) -> Dict:
    # if the host is not yet saved
    if host._state.adding:
        if hasattr(host, "_feature_sets"):
            return host._feature_sets
        else:
            return {}
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


def print_features(self: Data) -> str:
    from .._from_values import _print_values

    msg = ""
    features_lookup = Feature.lookup()
    for slot, feature_set in self.features._feature_set_by_slot.items():
        if feature_set.registry != "core.Feature":
            key_split = feature_set.registry.split(".")
            orm_name_with_schema = f"{key_split[0]}.{key_split[1]}"
            feature_sets_related_models = dict_related_model_to_related_name(
                feature_set
            )
            related_name = feature_sets_related_models.get(orm_name_with_schema)
            # first 5 feature records
            features = feature_set.__getattribute__(related_name).all()
            name_field = get_default_str_field(features[0])
            feature_names = [getattr(feature, name_field) for feature in features]
            host_id_field = get_host_id_field(self)
            kwargs = {host_id_field: self.id, "feature_set_id": feature_set.id}
            slots = self.feature_sets.through.objects.filter(**kwargs).list("slot")
            for slot in slots:
                msg += f"  {colors.bold(slot)}: {feature_set}\n"
                print_values = _print_values(feature_names, n=10)
                msg += f"    {print_values}\n"
        else:
            df_slot = feature_set.features.df()
            msg += f"  {colors.bold(slot)}: {feature_set}\n"
            for _, row in df_slot.iterrows():
                if row["type"] == "category" and row["registries"] is not None:
                    labels = self.labels.get(
                        getattr(features_lookup, row["name"]), mute=True
                    )
                    indent = ""
                    if isinstance(labels, dict):
                        msg += f"    🔗 {row['name']} ({row.registries})\n"
                        indent = "    "
                    else:
                        labels = {row["registries"]: labels}
                    for registry, labels in labels.items():
                        count_str = f"{len(labels)}, {colors.italic(f'{registry}')}"
                        field = get_default_str_field(labels)
                        print_values = _print_values(labels.list(field), n=10)
                        msg_objects = (
                            f"{indent}    🔗 {row['name']} ({count_str}):"
                            f" {print_values}\n"
                        )
                        msg += msg_objects
                else:
                    msg += f"    {row['name']} ({row['type']})\n"
    if msg != "":
        msg = f"{colors.green('Features')}:\n" + msg
    return msg


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
            return print_features(self._host)
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
