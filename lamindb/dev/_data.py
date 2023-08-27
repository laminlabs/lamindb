from collections import defaultdict
from typing import Any, Dict, List, Optional, Union

from lamin_utils import colors, logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Run
from lnschema_core.models import (
    Data,
    Dataset,
    Feature,
    FeatureSet,
    File,
    Label,
    Registry,
    format_field_value,
)

from .._query_set import QuerySet
from .._registry import get_default_str_field
from ._feature_manager import (
    FeatureManager,
    create_features_df,
    get_feature_set_links,
    get_host_id_field,
    get_label_links,
)
from ._run_context import run_context
from ._settings import settings
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


def describe(self):
    model_name = colors.green(self.__class__.__name__)
    msg = ""

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
    }
    if len(foreign_key_fields) > 0:
        record_msg = f"{model_name}({''.join([f'{i}={format_field_value(self.__getattribute__(i))}, ' for i in direct_fields])})"  # noqa
        msg += f"{record_msg.rstrip(', )')})\n\n"

        msg += f"{colors.green('Provenance')}:\n    "
        related_msg = "".join(
            [
                f"{emojis.get(i, '📎')} {i}: {self.__getattribute__(i)}\n    "
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

    feature_sets_all = self.feature_sets.all()
    if len(feature_sets_all) > 0:
        msg += f"{colors.green('Features')}:\n"
    feature_sets_subset = feature_sets_all.exclude(registry="core.Feature")
    if feature_sets_subset.exists():
        feature_sets_related_models = dict_related_model_to_related_name(
            feature_sets_subset.first()
        )
        for feature_set in feature_sets_subset:
            key_split = feature_set.registry.split(".")
            if len(key_split) == 3:
                logger.warning(
                    "you have a legacy entry in feature_set.field, should be just"
                    " 'bionty.Gene'"
                )
            orm_name_with_schema = f"{key_split[0]}.{key_split[1]}"
            field_name = "id"
            related_name = feature_sets_related_models.get(orm_name_with_schema)
            values = feature_set.__getattribute__(related_name).all()[:5].list("id")
            host_id_field = get_host_id_field(self)
            kwargs = {host_id_field: self.id, "feature_set_id": feature_set.id}
            slots = self.feature_sets.through.objects.filter(**kwargs).list("slot")
            for slot in slots:
                if slot == "var":
                    slot += " (X)"
                msg += f"  {colors.bold(slot)}:\n"
                ref = colors.italic(f"{orm_name_with_schema}.{field_name}")
                msg += f"    🔗 index ({feature_set.n}, {ref}): {values}\n".replace(
                    "]", "...]"
                )

    # display core.Feature features
    feature_sets = feature_sets_all.filter(registry="core.Feature").all()
    if feature_sets.exists():
        features_df = create_features_df(
            host=self, feature_sets=feature_sets.all(), exclude=True
        )
        for slot in features_df["slot"].unique():
            df_slot = features_df[features_df.slot == slot]
            if slot == "obs":
                slot += " (metadata)"
            msg += f"  {colors.bold(slot)}:\n"
            for _, row in df_slot.iterrows():
                labels = self.get_labels(row["name"], mute=True)
                indent = ""
                if isinstance(labels, dict):
                    msg += f"    🔗 {row['name']} ({row.registries})\n"
                    indent = "    "
                else:
                    labels = {row["registries"]: labels}
                for registry, labels in labels.items():
                    count_str = f"{len(labels)}, {colors.italic(f'{registry}')}"
                    try:
                        field = get_default_str_field(labels)
                    except ValueError:
                        field = "id"
                    values = labels.list(field)[:5]
                    msg_objects = (
                        f"{indent}    🔗 {row['name']} ({count_str}): {values}\n"
                    )
                    msg += msg_objects
    verbosity = settings.verbosity
    settings.verbosity = 3
    logger.info(msg)
    settings.verbosity = verbosity


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
    feature: Union[str, Registry],
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
            links_to_labels = get_label_links(self, registry, feature)
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
        raise ValueError("Please save the file/dataset before adding a label!")
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
            raise ValueError(
                "Please pass feature: add_labels(labels, feature='myfeature')"
            )
        else:
            record_feature = feature
        records_by_feature_orm[
            (record_feature, record.__class__.__get_name_with_schema__())
        ].append(record)
    for (feature, orm_name), records in records_by_feature_orm.items():
        getattr(self, self.features._accessor_by_orm[orm_name]).add(
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
            if "external" in linked_features_by_slot:
                feature_set = self.features._feature_set_by_slot["external"]
                features_list = feature_set.features.list()
            else:
                features_list = []
            features_list.append(feature)
            feature_set = FeatureSet(features_list, modality="meta")
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


@property  # type: ignore
@doc_args(Data.features.__doc__)
def features(self) -> "FeatureManager":
    """{}"""
    from lamindb.dev._feature_manager import FeatureManager

    return FeatureManager(self)


setattr(Data, "features", features)
setattr(Data, "add_labels", add_labels)
setattr(Data, "get_labels", get_labels)
setattr(Data, "describe", describe)
