from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Any, Iterable, List

from lamin_utils import colors, logger
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import (
    Artifact,
    Collection,
    Data,
    Feature,
    FeatureSet,
    Registry,
    Run,
    ULabel,
    __repr__,
    format_field_value,
)

from lamindb._feature_set import (
    dict_related_model_to_related_name,
    dict_schema_name_to_model_name,
)
from lamindb._parents import view_lineage
from lamindb._query_set import QuerySet
from lamindb.core._settings import settings

from ._feature_manager import (
    FeatureManager,
    get_feature_set_links,
    get_host_id_field,
    get_label_links,
    print_features,
)
from ._label_manager import LabelManager, print_labels
from ._run_context import run_context
from .exceptions import ValidationError

if TYPE_CHECKING:
    from lnschema_core.types import StrField

WARNING_RUN_TRANSFORM = "no run & transform get linked, consider calling ln.track()"


def get_run(run: Run | None) -> Run | None:
    if run is None:
        run = run_context.run
        if run is None and not settings.silence_file_run_transform_warning:
            logger.warning(WARNING_RUN_TRANSFORM)
    # suppress run by passing False
    elif not run:
        run = None
    return run


def add_transform_to_kwargs(kwargs: dict[str, Any], run: Run):
    if run is not None:
        kwargs["transform"] = run.transform


def save_feature_sets(self: Artifact | Collection) -> None:
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


def save_feature_set_links(self: Artifact | Collection) -> None:
    from lamindb._save import bulk_create

    Data = self.__class__
    if hasattr(self, "_feature_sets"):
        links = []
        host_id_field = get_host_id_field(self)
        for slot, feature_set in self._feature_sets.items():
            kwargs = {
                host_id_field: self.id,
                "feature_set_id": feature_set.id,
                "slot": slot,
            }
            links.append(Data.feature_sets.through(**kwargs))
        bulk_create(links, ignore_conflicts=True)


def format_repr(value: Registry, exclude: list[str] | str | None = None) -> str:
    if isinstance(exclude, str):
        exclude = [exclude]
    exclude_fields = set() if exclude is None else set(exclude)
    exclude_fields.update(["created_at", "updated_at"])

    fields = [
        f
        for f in value.__repr__(include_foreign_keys=False).split(", ")
        if not any(f"{excluded_field}=" in f for excluded_field in exclude_fields)
    ]
    repr = ", ".join(fields)
    if not repr.endswith(")"):
        repr += ")"
    return repr


@doc_args(Data.describe.__doc__)
def describe(self: Data):
    """{}."""
    # prefetch all many-to-many relationships
    # doesn't work for describing using artifact
    # self = (
    #     self.__class__.objects.using(self._state.db)
    #     .prefetch_related(
    #         *[f.name for f in self.__class__._meta.get_fields() if f.many_to_many]
    #     )
    #     .get(id=self.id)
    # )

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
    if not self._state.adding:
        # prefetch foreign key relationships
        self = (
            self.__class__.objects.using(self._state.db)
            .select_related(*foreign_key_fields)
            .get(id=self.id)
        )
        # prefetch m-2-m relationships
        self = (
            self.__class__.objects.using(self._state.db)
            .prefetch_related("feature_sets", "input_of")
            .get(id=self.id)
        )

    # provenance
    if len(foreign_key_fields) > 0:  # always True for Artifact and Collection
        record_msg = f"{colors.green(model_name)}{__repr__(self, include_foreign_keys=False).lstrip(model_name)}"
        msg += f"{record_msg}\n\n"

        msg += f"{colors.green('Provenance')}:\n  "
        related_msg = "".join(
            [
                f"📎 {field}: {format_repr(self.__getattribute__(field))}\n  "
                for field in foreign_key_fields
                if self.__getattribute__(field) is not None
            ]
        )
        msg += related_msg
    # input of
    if self.id is not None and self.input_of.exists():
        values = [format_field_value(i.started_at) for i in self.input_of.all()]
        msg += f"📎 input_of ({colors.italic('core.Run')}): {values}\n    "
    msg = msg.rstrip(" ")  # do not use removesuffix as we need to remove 2 or 4 spaces
    msg += print_features(self)
    msg += print_labels(self)

    logger.print(msg)


def validate_feature(feature: Feature, records: list[Registry]) -> None:
    """Validate feature record, adjust feature.type based on labels records."""
    if not isinstance(feature, Feature):
        raise TypeError("feature has to be of type Feature")
    if feature._state.adding:
        registries = {record.__class__.__get_name_with_schema__() for record in records}
        registries_str = "|".join(registries)
        msg = f"ln.Feature(name='{feature.name}', type='cat[{registries_str}]').save()"
        raise ValidationError(f"Feature not validated. If it looks correct: {msg}")


def get_labels(
    self,
    feature: Feature,
    mute: bool = False,
    flat_names: bool = False,
) -> QuerySet | dict[str, QuerySet] | list:
    """{}."""
    if not isinstance(feature, Feature):
        raise TypeError("feature has to be of type Feature")
    if feature.type is None or not feature.type.startswith("cat["):
        raise ValueError("feature does not have linked labels")
    registries_to_check = feature.type.replace("cat[", "").rstrip("]").split("|")
    if len(registries_to_check) > 1 and not mute:
        logger.warning("labels come from multiple registries!")
    # return an empty query set if self.id is still None
    if self.id is None:
        return QuerySet(self.__class__)
    qs_by_registry = {}
    for registry in registries_to_check:
        # currently need to distinguish between ULabel and non-ULabel, because
        # we only have the feature information for Label
        if registry == "core.ULabel":
            links_to_labels = get_label_links(self, registry, feature)
            label_ids = [link.ulabel_id for link in links_to_labels]
            qs_by_registry[registry] = ULabel.objects.using(self._state.db).filter(
                id__in=label_ids
            )
        else:
            qs_by_registry[registry] = getattr(
                self, self.features.accessor_by_orm[registry]
            ).all()
    if flat_names:
        # returns a flat list of names
        from lamindb._registry import get_default_str_field

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
    records: Registry | list[Registry] | QuerySet | Iterable,
    feature: Feature | None = None,
    *,
    field: StrField | None = None,
) -> None:
    """{}."""
    if self._state.adding:
        raise ValueError("Please save the artifact/collection before adding a label!")

    if isinstance(records, (QuerySet, QuerySet.__base__)):  # need to have both
        records = records.list()
    if isinstance(records, (str, Registry)):
        records = [records]
    if not isinstance(records, List):  # avoids warning for pd Series
        records = list(records)
    # create records from values
    if len(records) == 0:
        return None
    if isinstance(records[0], str):  # type: ignore
        records_validated = []
        # feature is needed if we want to create records from values
        if feature is None:
            raise ValueError(
                "Please pass a feature, e.g., via: label = ln.ULabel(name='my_label',"
                " feature=ln.Feature(name='my_feature'))"
            )
        if feature.type.startswith("cat["):
            orm_dict = dict_schema_name_to_model_name(Artifact)
            for reg in feature.type.replace("cat[", "").rstrip("]").split("|"):
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

    for record in records:
        if record._state.adding:
            raise ValidationError(
                f"{record} not validated. If it looks correct: record.save()"
            )

    if feature is None:
        d = dict_related_model_to_related_name(self.__class__)
        # strategy: group records by registry to reduce number of transactions
        records_by_related_name: dict = {}
        for record in records:
            related_name = d.get(record.__class__.__get_name_with_schema__())
            if related_name is None:
                raise ValueError(f"Can't add labels to {record.__class__} record!")
            if related_name not in records_by_related_name:
                records_by_related_name[related_name] = []
            records_by_related_name[related_name].append(record)
        for related_name, records in records_by_related_name.items():
            getattr(self, related_name).add(*records)
    else:
        validate_feature(feature, records)  # type:ignore
        records_by_registry = defaultdict(list)
        for record in records:
            records_by_registry[record.__class__.__get_name_with_schema__()].append(
                record
            )
        for registry_name, records in records_by_registry.items():
            labels_accessor = getattr(
                self, self.features.accessor_by_orm[registry_name]
            )
            # remove labels that are already linked as add doesn't perform update
            linked_labels = [r for r in records if r in labels_accessor.filter()]
            if len(linked_labels) > 0:
                labels_accessor.remove(*linked_labels)
            labels_accessor.add(*records, through_defaults={"feature_id": feature.id})
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
        for registry_name, _ in records_by_registry.items():
            msg = ""
            if not feature.type.startswith("cat[") or registry_name not in feature.type:
                if len(msg) > 0:
                    msg += ", "
                msg += f"linked feature '{feature.name}' to registry '{registry_name}'"
                if not feature.type.startswith("cat["):
                    feature.type = f"cat[{registry_name}]"
                elif registry_name not in feature.type:
                    feature.type = feature.type.rstrip("]") + f"|{registry_name}]"
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
                    feature_set = self.features.feature_set_by_slot["external"]
                    features_list = feature_set.features.list()
                else:
                    features_list = []
                features_list.append(feature)
                feature_set = FeatureSet(features_list)
                feature_set.save()
                if "external" in linked_features_by_slot:
                    old_feature_set_link = feature_set_links.filter(
                        slot="external"
                    ).one()
                    old_feature_set_link.delete()
                    remaining_links = self.feature_sets.through.objects.filter(
                        feature_set_id=feature_set.id
                    ).all()
                    if len(remaining_links) == 0:
                        old_feature_set = FeatureSet.filter(
                            id=old_feature_set_link.feature_set_id
                        ).one()
                        logger.info(
                            "nothing links to it anymore, deleting feature set"
                            f" {old_feature_set}"
                        )
                        old_feature_set.delete()
                self.features.add_feature_set(feature_set, slot="external")
                logger.save(
                    f"linked new feature '{feature.name}' together with new feature set"
                    f" {feature_set}"
                )


def _track_run_input(
    data: Data | Iterable[Data],
    is_run_input: bool | None = None,
    run: Run | None = None,
):
    # this is an internal hack right now for project-flow, but we can allow this
    # for the user in the future
    if isinstance(is_run_input, Run):
        run = is_run_input
        is_run_input = True
    elif run is None:
        run = run_context.run
    # consider that data is an iterable of Data
    data_iter: Iterable[Data] = [data] if isinstance(data, Data) else data
    track_run_input = False
    input_data = []
    if run is not None:
        # avoid cycles: data can't be both input and output
        def is_valid_input(data: Data):
            return (
                data.run_id != run.id
                and not data._state.adding
                and data._state.db in {"default", None}
            )

        input_data = [data for data in data_iter if is_valid_input(data)]
        input_data_ids = [data.id for data in input_data]
    if input_data:
        data_class_name = input_data[0].__class__.__name__.lower()
    # let us first look at the case in which the user does not
    # provide a boolean value for `is_run_input`
    # hence, we need to determine whether we actually want to
    # track a run or not
    if is_run_input is None:
        # we don't have a run record
        if run is None:
            if settings.track_run_inputs:
                logger.hint(
                    "you can auto-track these data as a run input by calling"
                    " `ln.track()`"
                )
        # assume we have a run record
        else:
            # assume there is non-cyclic candidate input data
            if input_data:
                if settings.track_run_inputs:
                    transform_note = ""
                    if len(input_data) == 1:
                        if input_data[0].transform is not None:
                            transform_note = (
                                ", adding parent transform"
                                f" {input_data[0].transform.id}"
                            )
                    logger.info(
                        f"adding {data_class_name} {input_data_ids} as input for run"
                        f" {run.id}{transform_note}"
                    )
                    track_run_input = True
                else:
                    logger.hint(
                        "track these data as a run input by passing `is_run_input=True`"
                    )
    else:
        track_run_input = is_run_input
    if track_run_input:
        if run is None:
            raise ValueError(
                "No run context set. Call ln.track() or link input to a"
                " run object via `run.input_artifacts.add(artifact)`"
            )
        # avoid adding the same run twice
        run.save()
        if data_class_name == "artifact":
            LinkORM = run.input_artifacts.through
            links = [
                LinkORM(run_id=run.id, artifact_id=data_id)
                for data_id in input_data_ids
            ]
        else:
            LinkORM = run.input_collections.through
            links = [
                LinkORM(run_id=run.id, collection_id=data_id)
                for data_id in input_data_ids
            ]
        LinkORM.objects.bulk_create(links, ignore_conflicts=True)
        # generalize below for more than one data batch
        if len(input_data) == 1:
            if input_data[0].transform is not None:
                run.transform.parents.add(input_data[0].transform)


@property  # type: ignore
@doc_args(Data.features.__doc__)
def features(self) -> FeatureManager:
    """{}."""
    from lamindb.core._feature_manager import FeatureManager

    return FeatureManager(self)


@property  # type: ignore
@doc_args(Data.labels.__doc__)
def labels(self) -> LabelManager:
    """{}."""
    from lamindb.core._label_manager import LabelManager

    return LabelManager(self)


Data.features = features
Data.labels = labels
Data.describe = describe
Data.view_lineage = view_lineage
