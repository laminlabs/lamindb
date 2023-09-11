from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Union

from lamin_utils import colors, logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core.models import (
    Data,
    Dataset,
    Feature,
    FeatureSet,
    File,
    Registry,
    Run,
    ULabel,
    __repr__,
    format_field_value,
)
from lnschema_core.types import StrField

from .._parents import view_flow_dataset, view_flow_file
from .._query_set import QuerySet
from ._feature_manager import (
    FeatureManager,
    dict_schema_name_to_model_name,
    get_feature_set_links,
    get_host_id_field,
    get_label_links,
    print_features,
)
from ._label_manager import LabelManager, print_labels
from ._priors import priors
from ._run_context import run_context
from .exceptions import ValidationError


def get_run(run: Optional[Run]) -> Optional[Run]:
    if run is None:
        run = run_context.run
        if run is None:
            logger.warning(
                "no run & transform get linked, consider passing a `run` or calling"
                " ln.track()"
            )
    return run


def add_transform_to_kwargs(kwargs: Dict[str, Any], run: Run):
    if run is not None:
        kwargs["transform"] = run.transform


def save_transform_run_feature_sets(self: Union[File, Dataset]) -> None:
    if self.transform is not None:
        self.transform.save()
    if self.run is not None:
        self.run.save()
    if hasattr(self, "_feature_sets"):
        saved_feature_sets = {}
        for key, feature_set in self._feature_sets.items():
            if isinstance(feature_set, FeatureSet) and feature_set._state.adding:
                feature_set.save()
                saved_feature_sets[key] = feature_set
        if len(saved_feature_sets) > 0:
            s = "s" if len(saved_feature_sets) > 1 else ""
            display_feature_set_keys = ",".join(
                f"'{key}'" for key in saved_feature_sets.keys()
            )
            logger.save(
                f"saved {len(saved_feature_sets)} feature set{s} for slot{s}:"
                f" {display_feature_set_keys}"
            )


def save_feature_set_links(self: Union[File, Dataset]) -> None:
    from lamindb._save import bulk_create

    Data = self.__class__
    if hasattr(self, "_feature_sets"):
        links = []
        host_id_field = get_host_id_field(self)
        for slot, feature_set in self._feature_sets.items():
            if isinstance(feature_set, str):
                assert len(feature_set) == 20
                feature_set_id = feature_set
            else:
                feature_set_id = feature_set.id
            kwargs = {
                host_id_field: self.id,
                "feature_set_id": feature_set_id,
                "slot": slot,
            }
            links.append(Data.feature_sets.through(**kwargs))
        bulk_create(links)


@doc_args(Data.describe.__doc__)
def describe(self: Data):
    """{}"""
    model_name = self.__class__.__name__
    msg = ""

    fields = self._meta.fields
    direct_fields = []
    foreign_key_fields = []
    for f in fields:
        if f.is_relation:
            foreign_key_fields.append(f.name)
        else:
            direct_fields.append(f.name)

    # Display Provenance
    # display line by line the foreign key fields
    from .._parents import _transform_emoji

    emojis = {
        "storage": "🗃️",
        "created_by": "👤",
        "transform": _transform_emoji(self.transform),
        "run": "👣",
        "initial_version": "🔖",
        "file": "📁",
    }
    if len(foreign_key_fields) > 0:  # always True for File and Dataset
        record_msg = f"{colors.green(model_name)}{__repr__(self, include_foreign_keys=False).lstrip(model_name)}"  # noqa
        msg += f"{record_msg}\n\n"

        msg += f"{colors.green('Provenance')}:\n  "
        related_msg = "".join(
            [
                f"{emojis.get(i, '📎')} {i}: {self.__getattribute__(i)}\n  "
                for i in foreign_key_fields
                if self.__getattribute__(i) is not None
            ]
        )
        msg += related_msg
    # input of
    if self.input_of.exists():
        values = [format_field_value(i.run_at) for i in self.input_of.all()]
        msg += f"⬇️ input_of ({colors.italic('core.Run')}): {values}\n    "
    msg = msg.rstrip("    ")
    msg += print_features(self)
    msg += print_labels(self)

    logger.print(msg)


def validate_feature(feature: Feature, records: List[Registry]) -> None:
    if not isinstance(feature, Feature):
        raise TypeError("feature has to be of type Feature")
    if feature._state.adding:
        registries = set(
            [record.__class__.__get_name_with_schema__() for record in records]
        )
        registries_str = "|".join(registries)
        msg = (
            f"ln.Feature(name='{feature.name}', type='category',"
            f" registries='{registries_str}').save()"
        )
        raise ValidationError(f"Feature not validated. If it looks correct: {msg}")


def get_labels(
    self,
    feature: Feature,
    mute: bool = False,
    flat_names: bool = False,
) -> Union[QuerySet, Dict[str, QuerySet], List]:
    """{}"""
    if not isinstance(feature, Feature):
        raise TypeError("feature has to be of type Feature")
    if feature.registries is None:
        raise ValueError("feature does not have linked labels")
    registries_to_check = feature.registries.split("|")
    if len(registries_to_check) > 1 and not mute:
        logger.warning("labels come from multiple registries!")
    qs_by_registry = {}
    for registry in registries_to_check:
        # currently need to distinguish between Label and non-Label, because
        # we only have the feature information for Label
        if registry == "core.ULabel":
            links_to_labels = get_label_links(self, registry, feature)
            label_ids = [link.ulabel_id for link in links_to_labels]
            qs_by_registry[registry] = ULabel.objects.filter(id__in=label_ids)
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


def add_labels(
    self,
    records: Union[Registry, List[Registry], QuerySet, Iterable],
    feature: Feature,
    *,
    field: Optional[StrField] = None,
) -> None:
    """{}"""
    if isinstance(records, (QuerySet, QuerySet.__base__)):  # need to have both
        records = records.list()
    if isinstance(records, (str, Registry)):
        records = [records]
    if not isinstance(records, List):  # avoids warning for pd Series
        records = list(records)
    if isinstance(records[0], str):  # type: ignore
        records_validated = []
        if feature.registries is not None:
            orm_dict = dict_schema_name_to_model_name(File)
            for reg in feature.registries.split("|"):
                orm = orm_dict.get(reg)
                records_validated += orm.from_values(records, field=field)

        # feature doesn't have registries and therefore can't create records from values
        # ask users to pass records
        if len(records_validated) == 0:
            raise ValueError(
                "Please pass a record (a `Registry` object), not a string, e.g., via:"
                " label"
                f" = ln.ULabel(name='{records[0]}')"  # type: ignore
            )
        records = records_validated

    if self._state.adding:
        raise ValueError("Please save the file/dataset before adding a label!")
    for record in records:
        if record._state.adding:
            raise ValidationError(
                f"{record} not validated. If it looks correct: record.save()"
            )
    validate_feature(feature, records)  # type:ignore
    records_by_registry = defaultdict(list)
    for record in records:
        records_by_registry[record.__class__.__get_name_with_schema__()].append(record)
    for registry_name, records in records_by_registry.items():
        getattr(self, self.features._accessor_by_orm[registry_name]).add(
            *records, through_defaults={"feature_id": feature.id}
        )
    feature_set_links = get_feature_set_links(self)
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
    for registry_name, records in records_by_registry.items():
        msg = ""
        if feature.registries is None or registry_name not in feature.registries:
            if len(msg) > 0:
                msg += ", "
            msg += f"linked feature '{feature.name}' to registry '{registry_name}'"
            if feature.registries is None:
                feature.registries = registry_name
            elif registry_name not in feature.registries:
                feature.registries += f"|{registry_name}"
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
            if "external" in linked_features_by_slot:
                feature_set = self.features._feature_set_by_slot["external"]
                features_list = feature_set.features.list()
            else:
                features_list = []
            features_list.append(feature)
            feature_set = FeatureSet(features_list, modality=priors.modalities.meta)
            feature_set.save()
            if "external" in linked_features_by_slot:
                old_feature_set_link = feature_set_links.filter(slot="external").one()
                old_feature_set_link.delete()
                remaining_links = self.feature_sets.through.objects.filter(
                    feature_set_id=feature_set.id
                ).all()
                if len(remaining_links) == 0:
                    old_feature_set = FeatureSet.filter(
                        id=old_feature_set_link.feature_set_id
                    ).one()
                    logger.info(
                        "no file links to it anymore, deleting feature set"
                        f" {old_feature_set}"
                    )
                    old_feature_set.delete()
            self.features.add_feature_set(feature_set, slot="external")
            logger.save(
                f"linked new feature '{feature.name}' together with new feature set"
                f" {feature_set}"
            )


def view_flow(data: Data, with_children: bool = True) -> None:
    if isinstance(data, File):
        view_flow_file(file=data, with_children=with_children)
    elif isinstance(data, Dataset):
        view_flow_dataset(dataset=data, with_children=with_children)


@property  # type: ignore
@doc_args(Data.features.__doc__)
def features(self) -> "FeatureManager":
    """{}"""
    from lamindb.dev._feature_manager import FeatureManager

    return FeatureManager(self)


@property  # type: ignore
@doc_args(Data.labels.__doc__)
def labels(self) -> "LabelManager":
    """{}"""
    from lamindb.dev._label_manager import LabelManager

    return LabelManager(self)


setattr(Data, "features", features)
setattr(Data, "labels", labels)
setattr(Data, "describe", describe)
setattr(Data, "view_flow", view_flow)
