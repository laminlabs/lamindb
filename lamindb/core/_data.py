from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Any

from django.db import connections
from lamin_utils import colors, logger
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import (
    Artifact,
    Collection,
    Feature,
    FeatureSet,
    Record,
    Run,
    ULabel,
    format_field_value,
    record_repr,
)

from lamindb._parents import view_lineage
from lamindb._query_set import QuerySet
from lamindb._record import get_name_field
from lamindb.core._settings import settings

from ._context import context
from ._django import get_artifact_with_related, get_related_model
from ._feature_manager import (
    get_feature_set_links,
    get_host_id_field,
    get_label_links,
    print_features,
)
from ._label_manager import print_labels
from .exceptions import ValidationError
from .schema import (
    dict_related_model_to_related_name,
    dict_schema_name_to_model_name,
)

if TYPE_CHECKING:
    from collections.abc import Iterable

    from lnschema_core.types import StrField


WARNING_RUN_TRANSFORM = "no run & transform got linked, call `ln.track()` & re-run"

WARNING_NO_INPUT = "run input wasn't tracked, call `ln.track()` and re-run"


def get_run(run: Run | None) -> Run | None:
    if run is None:
        run = context.run
        if run is None and not settings.creation.artifact_silence_missing_run_warning:
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
        from lamindb.core._feature_manager import get_feature_set_by_slot_

        existing_feature_sets = get_feature_set_by_slot_(self)
        saved_feature_sets = {}
        for key, feature_set in self._feature_sets.items():
            if isinstance(feature_set, FeatureSet) and feature_set._state.adding:
                feature_set.save()
                saved_feature_sets[key] = feature_set
            if key in existing_feature_sets:
                # remove existing feature set on the same slot
                self.feature_sets.remove(existing_feature_sets[key])
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
                "featureset_id": feature_set.id,
                "slot": slot,
            }
            links.append(Data.feature_sets.through(**kwargs))
        bulk_create(links, ignore_conflicts=True)


def format_provenance(self, fk_data, print_types):
    type_str = lambda attr: (
        f": {get_related_model(self.__class__, attr).__name__}" if print_types else ""
    )

    return "".join(
        [
            f"    .{field_name}{type_str(field_name)} = {format_field_value(value.get('name'))}\n"
            for field_name, value in fk_data.items()
            if value.get("name")
        ]
    )


def format_input_of_runs(self, print_types):
    if self.id is not None and self.input_of_runs.exists():
        values = [format_field_value(i.started_at) for i in self.input_of_runs.all()]
        type_str = ": Run" if print_types else ""  # type: ignore
        return f"    .input_of_runs{type_str} = {', '.join(values)}\n"
    return ""


def format_labels_and_features(self, related_data, print_types):
    msg = print_labels(
        self, m2m_data=related_data.get("m2m", {}), print_types=print_types
    )
    if isinstance(self, Artifact):
        msg += print_features(  # type: ignore
            self,
            related_data=related_data,
            print_types=print_types,
            print_params=hasattr(self, "type") and self.type == "model",
        )
    return msg


def _describe_postgres(self: Artifact | Collection, print_types: bool = False):
    model_name = self.__class__.__name__
    msg = f"{colors.green(model_name)}{record_repr(self, include_foreign_keys=False).lstrip(model_name)}\n"
    if self._state.db is not None and self._state.db != "default":
        msg += f"  {colors.italic('Database instance')}\n"
        msg += f"    slug: {self._state.db}\n"

    if model_name == "Artifact":
        result = get_artifact_with_related(
            self,
            include_feature_link=True,
            include_fk=True,
            include_m2m=True,
            include_featureset=True,
        )
    else:
        result = get_artifact_with_related(self, include_fk=True, include_m2m=True)
    related_data = result.get("related_data", {})
    fk_data = related_data.get("fk", {})

    # Provenance
    prov_msg = format_provenance(self, fk_data, print_types)
    if prov_msg:
        msg += f"  {colors.italic('Provenance')}\n{prov_msg}"

    # Input of runs
    input_of_message = format_input_of_runs(self, print_types)
    if input_of_message:
        msg += f"  {colors.italic('Usage')}\n{input_of_message}"

    # Labels and features
    msg += format_labels_and_features(self, related_data, print_types)

    # Print entire message
    logger.print(msg)


@doc_args(Artifact.describe.__doc__)
def describe(self: Artifact | Collection, print_types: bool = False):
    """{}"""  # noqa: D415
    if not self._state.adding and connections[self._state.db].vendor == "postgresql":
        return _describe_postgres(self, print_types=print_types)

    model_name = self.__class__.__name__
    msg = f"{colors.green(model_name)}{record_repr(self, include_foreign_keys=False).lstrip(model_name)}\n"
    if self._state.db is not None and self._state.db != "default":
        msg += f"  {colors.italic('Database instance')}\n"
        msg += f"    slug: {self._state.db}\n"

    prov_msg = ""
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
        many_to_many_fields = []
        if isinstance(self, (Collection, Artifact)):
            many_to_many_fields.append("input_of_runs")
        if isinstance(self, Artifact):
            many_to_many_fields.append("feature_sets")
        self = (
            self.__class__.objects.using(self._state.db)
            .prefetch_related(*many_to_many_fields)
            .get(id=self.id)
        )

    # provenance
    if len(foreign_key_fields) > 0:  # always True for Artifact and Collection
        fields_values = [(field, getattr(self, field)) for field in foreign_key_fields]
        type_str = lambda attr: (
            f": {attr.__class__.__get_name_with_schema__()}" if print_types else ""
        )
        related_msg = "".join(
            [
                f"    .{field_name}{type_str(attr)} = {format_field_value(getattr(attr, get_name_field(attr)))}\n"
                for (field_name, attr) in fields_values
                if attr is not None
            ]
        )
        prov_msg += related_msg
    if prov_msg:
        msg += f"  {colors.italic('Provenance')}\n"
        msg += prov_msg

    # Input of runs
    input_of_message = format_input_of_runs(self, print_types)
    if input_of_message:
        msg += f"  {colors.italic('Usage')}\n{input_of_message}"

    # Labels and features
    msg += format_labels_and_features(self, {}, print_types)

    # Print entire message
    logger.print(msg)


def validate_feature(feature: Feature, records: list[Record]) -> None:
    """Validate feature record, adjust feature.dtype based on labels records."""
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
    """{}"""  # noqa: D415
    if not isinstance(feature, Feature):
        raise TypeError("feature has to be of type Feature")
    if feature.dtype is None or not feature.dtype.startswith("cat["):
        raise ValueError("feature does not have linked labels")
    registries_to_check = feature.dtype.replace("cat[", "").rstrip("]").split("|")
    if len(registries_to_check) > 1 and not mute:
        logger.warning("labels come from multiple registries!")
    # return an empty query set if self.id is still None
    if self.id is None:
        return QuerySet(self.__class__)
    qs_by_registry = {}
    for registry in registries_to_check:
        # currently need to distinguish between ULabel and non-ULabel, because
        # we only have the feature information for Label
        if registry == "ULabel":
            links_to_labels = get_label_links(self, registry, feature)
            label_ids = [link.ulabel_id for link in links_to_labels]
            qs_by_registry[registry] = ULabel.objects.using(self._state.db).filter(
                id__in=label_ids
            )
        elif registry in self.features._accessor_by_registry:
            qs_by_registry[registry] = getattr(
                self, self.features._accessor_by_registry[registry]
            ).all()
    if flat_names:
        # returns a flat list of names
        from lamindb._record import get_name_field

        values = []
        for v in qs_by_registry.values():
            values += v.list(get_name_field(v))
        return values
    if len(registries_to_check) == 1 and registry in qs_by_registry:
        return qs_by_registry[registry]
    else:
        return qs_by_registry


def add_labels(
    self,
    records: Record | list[Record] | QuerySet | Iterable,
    feature: Feature | None = None,
    *,
    field: StrField | None = None,
    feature_ref_is_name: bool | None = None,
    label_ref_is_name: bool | None = None,
) -> None:
    """{}"""  # noqa: D415
    if self._state.adding:
        raise ValueError("Please save the artifact/collection before adding a label!")

    if isinstance(records, (QuerySet, QuerySet.__base__)):  # need to have both
        records = records.list()
    if isinstance(records, (str, Record)):
        records = [records]
    if not isinstance(records, list):  # avoids warning for pd Series
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
        if feature.dtype.startswith("cat["):
            orm_dict = dict_schema_name_to_model_name(Artifact)
            for reg in feature.dtype.replace("cat[", "").rstrip("]").split("|"):
                registry = orm_dict.get(reg)
                records_validated += registry.from_values(records, field=field)

        # feature doesn't have registries and therefore can't create records from values
        # ask users to pass records
        if len(records_validated) == 0:
            raise ValueError(
                "Please pass a record (a `Record` object), not a string, e.g., via:"
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
            if registry_name not in self.features._accessor_by_registry:
                logger.warning(f"skipping {registry_name}")
                continue
            if len(records) == 0:
                continue
            labels_accessor = getattr(
                self, self.features._accessor_by_registry[registry_name]
            )
            # remove labels that are already linked as add doesn't perform update
            linked_labels = [r for r in records if r in labels_accessor.filter()]
            if len(linked_labels) > 0:
                labels_accessor.remove(*linked_labels)
            through_defaults = {"feature_id": feature.id}
            link_model = labels_accessor.through
            if hasattr(link_model, "feature_ref_is_name"):
                through_defaults["feature_ref_is_name"] = feature_ref_is_name
                through_defaults["label_ref_is_name"] = label_ref_is_name
            labels_accessor.add(*records, through_defaults=through_defaults)
        links_feature_set = get_feature_set_links(self)
        feature_set_ids = [link.featureset_id for link in links_feature_set.all()]
        # get all linked features of type Feature
        feature_sets = FeatureSet.filter(id__in=feature_set_ids).all()
        {
            links_feature_set.filter(featureset_id=feature_set.id)
            .one()
            .slot: feature_set.features.all()
            for feature_set in feature_sets
            if "Feature" == feature_set.registry
        }
        for registry_name, _ in records_by_registry.items():
            if registry_name not in feature.dtype:
                logger.debug(
                    f"updated categorical feature '{feature.name}' type with registry '{registry_name}'"
                )
                if not feature.dtype.startswith("cat["):
                    feature.dtype = f"cat[{registry_name}]"
                elif registry_name not in feature.dtype:
                    feature.dtype = feature.dtype.rstrip("]") + f"|{registry_name}]"
                feature.save()


def _track_run_input(
    data: Artifact | Collection | Iterable[Artifact] | Iterable[Collection],
    is_run_input: bool | None = None,
    run: Run | None = None,
):
    # this is an internal hack right now for project-flow, but we can allow this
    # for the user in the future
    if isinstance(is_run_input, Run):
        run = is_run_input
        is_run_input = True
    elif run is None:
        run = context.run
    # consider that data is an iterable of Data
    data_iter: Iterable[Artifact] | Iterable[Collection] = (
        [data] if isinstance(data, (Artifact, Collection)) else data
    )
    track_run_input = False
    input_data = []
    if run is not None:
        # avoid cycles: data can't be both input and output
        def is_valid_input(data: Artifact | Collection):
            is_valid = False
            if data._state.db == "default":
                # things are OK if the record is on the default db
                is_valid = True
            elif data._state.db is None:
                # if a record is not yet saved, it can't be an input
                # we silently ignore because what likely happens is that
                # the user works with an object that's about to be saved
                # in the current Python session
                is_valid = False
            else:
                # record is on another db
                # we have to save the record into the current db with
                # the run being attached to a transfer transform
                logger.important(
                    f"completing transfer to track {data.__class__.__name__}('{data.uid[:8]}') as input"
                )
                data.save()
                is_valid = True
            return (
                data.run_id != run.id
                and not data._state.adding  # this seems duplicated with data._state.db is None
                and is_valid
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
                logger.warning(WARNING_NO_INPUT)
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
                        f"adding {data_class_name} ids {input_data_ids} as inputs for run"
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
            raise ValueError("No run context set. Call `ln.track()`.")
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
                run.transform.predecessors.add(input_data[0].transform)
