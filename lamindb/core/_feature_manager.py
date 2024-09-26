from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from itertools import compress
from typing import TYPE_CHECKING, Any

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData
from django.contrib.postgres.aggregates import ArrayAgg
from django.db import connections
from django.db.models import Aggregate
from lamin_utils import colors, logger
from lamindb_setup.core.upath import create_path
from lnschema_core.models import (
    Artifact,
    Collection,
    Feature,
    FeatureManager,
    FeatureValue,
    LinkORM,
    Param,
    ParamManager,
    ParamManagerArtifact,
    ParamManagerRun,
    ParamValue,
    Record,
    Run,
    ULabel,
)

from lamindb._feature import FEATURE_TYPES, convert_numpy_dtype_to_lamin_feature_type
from lamindb._feature_set import DICT_KEYS_TYPE, FeatureSet
from lamindb._record import (
    REGISTRY_UNIQUE_FIELD,
    get_name_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)
from lamindb._save import save
from lamindb.core.exceptions import ValidationError
from lamindb.core.storage import LocalPathClasses

from ._django import get_artifact_with_related
from ._settings import settings
from .schema import (
    dict_related_model_to_related_name,
)

if TYPE_CHECKING:
    from lnschema_core.types import FieldAttr

    from lamindb._query_set import QuerySet


def get_host_id_field(host: Artifact | Collection) -> str:
    if isinstance(host, Artifact):
        host_id_field = "artifact_id"
    else:
        host_id_field = "collection_id"
    return host_id_field


def get_accessor_by_registry_(host: Artifact | Collection) -> dict:
    dictionary = {
        field.related_model.__get_name_with_schema__(): field.name
        for field in host._meta.related_objects
    }
    dictionary["Feature"] = "features"
    dictionary["ULabel"] = "ulabels"
    return dictionary


def get_feature_set_by_slot_(host) -> dict:
    # if the host is not yet saved
    if host._state.adding:
        if hasattr(host, "_feature_sets"):
            return host._feature_sets
        else:
            return {}
    host_db = host._state.db
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id}
    # otherwise, we need a query
    links_feature_set = (
        host.feature_sets.through.objects.using(host_db)
        .filter(**kwargs)
        .select_related("featureset")
    )
    return {fsl.slot: fsl.featureset for fsl in links_feature_set}


def get_label_links(
    host: Artifact | Collection, registry: str, feature: Feature
) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id, "feature_id": feature.id}
    link_records = (
        getattr(host, host.features._accessor_by_registry[registry])
        .through.objects.using(host._state.db)
        .filter(**kwargs)
    )
    return link_records


def get_feature_set_links(host: Artifact | Collection) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id}
    links_feature_set = host.feature_sets.through.objects.filter(**kwargs)
    return links_feature_set


def get_link_attr(link: LinkORM | type[LinkORM], data: Artifact | Collection) -> str:
    link_model_name = link.__class__.__name__
    if link_model_name in {"Registry", "ModelBase"}:  # we passed the type of the link
        link_model_name = link.__name__
    return link_model_name.replace(data.__class__.__name__, "").lower()


# Custom aggregation for SQLite
class GroupConcat(Aggregate):
    function = "GROUP_CONCAT"
    template = '%(function)s(%(expressions)s, ", ")'


def custom_aggregate(field, using: str):
    if connections[using].vendor == "postgresql":
        return ArrayAgg(field)
    else:
        return GroupConcat(field)


def print_features(
    self: Artifact | Collection,
    related_data: dict | None = None,
    print_types: bool = False,
    to_dict: bool = False,
    print_params: bool = False,
) -> str | dict[str, Any]:
    from lamindb._from_values import _print_values

    if not related_data:
        artifact_meta = get_artifact_with_related(self, fks=False)
        related_data = artifact_meta.get("related_data", {})

    m2m_data = related_data.get("m2m", {}) if related_data else {}
    m2m_name = {}
    for related_name, values in m2m_data.items():
        link_model = getattr(self.__class__, related_name).through
        related_model_name = link_model.__name__.replace(
            self.__class__.__name__, ""
        ).lower()
        m2m_name[related_model_name] = values
    links_data = related_data.get("link", {}) if related_data else {}
    feature_dict = {
        id: (name, dtype)
        for id, name, dtype in Feature.objects.using(self._state.db).values_list(
            "id", "name", "dtype"
        )
    }

    msg = ""
    dictionary = {}

    # categorical feature values
    if not print_params:
        labels_msg = ""
        labels_msgs = []
        feature_values: dict = {}
        for link_name, link_values in links_data.items():
            related_name = link_name.removeprefix("links_").replace("_", "")
            link_model = getattr(self.__class__, link_name).rel.related_model
            if not link_values:
                continue
            for link_value in link_values:
                feature_id = link_value.get("feature")
                if feature_id is None:
                    continue
                feature_name = feature_dict.get(feature_id)[0]
                if feature_name not in feature_values:
                    feature_values[feature_name] = (feature_dict.get(feature_id)[1], [])
                label_id = link_value.get(related_name)
                feature_values[feature_name][1].append(
                    m2m_name.get(related_name, {}).get(label_id)
                )
        for feature_name, (dtype, labels_list) in feature_values.items():
            print_values = _print_values(labels_list, n=10)
            type_str = f": {dtype}" if print_types else ""
            if to_dict:
                dictionary[feature_name] = (
                    labels_list if len(labels_list) > 1 else labels_list[0]
                )
            labels_msgs.append(f"    '{feature_name}'{type_str} = {print_values}")
        if len(labels_msgs) > 0:
            labels_msg = "\n".join(sorted(labels_msgs)) + "\n"
            msg += labels_msg

    # non-categorical feature values
    non_labels_msg = ""
    if self.id is not None and self.__class__ == Artifact or self.__class__ == Run:
        attr_name = "param" if print_params else "feature"
        _feature_values = (
            getattr(self, f"_{attr_name}_values")
            .values(f"{attr_name}__name", f"{attr_name}__dtype")
            .annotate(values=custom_aggregate("value", self._state.db))
            .order_by(f"{attr_name}__name")
        )
        if len(_feature_values) > 0:
            for fv in _feature_values:
                feature_name = fv[f"{attr_name}__name"]
                feature_dtype = fv[f"{attr_name}__dtype"]
                values = fv["values"]
                # TODO: understand why the below is necessary
                if not isinstance(values, list):
                    values = [values]
                if to_dict:
                    dictionary[feature_name] = values if len(values) > 1 else values[0]
                type_str = f": {feature_dtype}" if print_types else ""
                printed_values = (
                    _print_values(values, n=10, quotes=False)
                    if not feature_dtype.startswith("list")
                    else values
                )
                non_labels_msg += f"    '{feature_name}'{type_str} = {printed_values}\n"
            msg += non_labels_msg

    if msg != "":
        header = "Features" if not print_params else "Params"
        msg = f"  {colors.italic(header)}\n" + msg

    # feature sets
    if not print_params:
        feature_set_msg = ""
        for slot, feature_set in get_feature_set_by_slot_(self).items():
            features = feature_set.members
            # features.first() is a lot slower than features[0] here
            name_field = get_name_field(features[0])
            feature_names = list(features.values_list(name_field, flat=True)[:20])
            type_str = f": {feature_set.registry}" if print_types else ""
            feature_set_msg += (
                f"    '{slot}'{type_str} = {_print_values(feature_names)}\n"
            )
        if feature_set_msg:
            msg += f"  {colors.italic('Feature sets')}\n"
            msg += feature_set_msg
    if to_dict:
        return dictionary
    else:
        return msg


def parse_feature_sets_from_anndata(
    adata: AnnData,
    var_field: FieldAttr | None = None,
    obs_field: FieldAttr = Feature.name,
    mute: bool = False,
    organism: str | Record | None = None,
) -> dict:
    data_parse = adata
    if not isinstance(adata, AnnData):  # is a path
        filepath = create_path(adata)  # returns Path for local
        if not isinstance(filepath, LocalPathClasses):
            from lamindb.core.storage._backed_access import backed_access

            using_key = settings._using_key
            data_parse = backed_access(filepath, using_key=using_key)
        else:
            data_parse = ad.read_h5ad(filepath, backed="r")
        type = "float"
    else:
        type = (
            "float"
            if adata.X is None
            else convert_numpy_dtype_to_lamin_feature_type(adata.X.dtype)
        )
    feature_sets = {}
    if var_field is not None:
        logger.info("parsing feature names of X stored in slot 'var'")
        logger.indent = "   "
        feature_set_var = FeatureSet.from_values(
            data_parse.var.index,
            var_field,
            type=type,
            mute=mute,
            organism=organism,
            raise_validation_error=False,
        )
        if feature_set_var is not None:
            feature_sets["var"] = feature_set_var
            logger.save(f"linked: {feature_set_var}")
        logger.indent = ""
        if feature_set_var is None:
            logger.warning("skip linking features to artifact in slot 'var'")
    if len(data_parse.obs.columns) > 0:
        logger.info("parsing feature names of slot 'obs'")
        logger.indent = "   "
        feature_set_obs = FeatureSet.from_df(
            df=data_parse.obs,
            field=obs_field,
            mute=mute,
            organism=organism,
        )
        if feature_set_obs is not None:
            feature_sets["obs"] = feature_set_obs
            logger.save(f"linked: {feature_set_obs}")
        logger.indent = ""
        if feature_set_obs is None:
            logger.warning("skip linking features to artifact in slot 'obs'")
    return feature_sets


def infer_feature_type_convert_json(
    value: Any, mute: bool = False, str_as_ulabel: bool = True
) -> tuple[str, Any]:
    if isinstance(value, bool):
        return FEATURE_TYPES["bool"], value
    elif isinstance(value, int):
        return FEATURE_TYPES["int"], value
    elif isinstance(value, float):
        return FEATURE_TYPES["float"], value
    elif isinstance(value, str):
        if str_as_ulabel:
            return FEATURE_TYPES["str"] + "[ULabel]", value
        else:
            return "str", value
    elif isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        if isinstance(value, (pd.Series, np.ndarray)):
            return convert_numpy_dtype_to_lamin_feature_type(
                value.dtype, str_as_cat=str_as_ulabel
            ), list(value)
        if isinstance(value, dict):
            return "dict", value
        if len(value) > 0:  # type: ignore
            first_element_type = type(next(iter(value)))
            if all(isinstance(elem, first_element_type) for elem in value):
                if first_element_type is bool:
                    return f"list[{FEATURE_TYPES['bool']}]", value
                elif first_element_type is int:
                    return f"list[{FEATURE_TYPES['int']}]", value
                elif first_element_type is float:
                    return f"list[{FEATURE_TYPES['float']}]", value
                elif first_element_type is str:
                    if str_as_ulabel:
                        return FEATURE_TYPES["str"] + "[ULabel]", value
                    else:
                        return "list[str]", value
                elif first_element_type == Record:
                    return (
                        f"cat[{first_element_type.__get_name_with_schema__()}]",
                        value,
                    )
    elif isinstance(value, Record):
        return (f"cat[{value.__class__.__get_name_with_schema__()}]", value)
    if not mute:
        logger.warning(f"cannot infer feature type of: {value}, returning '?")
    return ("?", value)


def __init__(self, host: Artifact | Collection | Run):
    self._host = host
    self._feature_set_by_slot_ = None
    self._accessor_by_registry_ = None


def __repr__(self) -> str:
    return print_features(self._host, print_params=(self.__class__ == ParamManager))  # type: ignore


def get_values(self) -> dict[str, Any]:
    """Get feature values as a dictionary."""
    return print_features(
        self._host, to_dict=True, print_params=(self.__class__ == ParamManager)
    )  # type: ignore


def __getitem__(self, slot) -> QuerySet:
    if slot not in self._feature_set_by_slot:
        raise ValueError(
            f"No linked feature set for slot: {slot}\nDid you get validation"
            " warnings? Only features that match registered features get validated"
            " and linked."
        )
    feature_set = self._feature_set_by_slot[slot]
    orm_name = feature_set.registry
    return getattr(feature_set, self._accessor_by_registry[orm_name]).all()


def filter_base(cls, **expression):
    if cls is FeatureManager:
        model = Feature
        value_model = FeatureValue
    else:
        model = Param
        value_model = ParamValue
    keys_normalized = [key.split("__")[0] for key in expression]
    validated = model.validate(keys_normalized, field="name", mute=True)
    if sum(validated) != len(keys_normalized):
        raise ValidationError(
            f"Some keys in the filter expression are not registered as features: {np.array(keys_normalized)[~validated]}"
        )
    new_expression = {}
    features = model.filter(name__in=keys_normalized).all().distinct()
    for key, value in expression.items():
        split_key = key.split("__")
        normalized_key = split_key[0]
        comparator = ""
        if len(split_key) == 2:
            comparator = f"__{split_key[1]}"
        feature = features.get(name=normalized_key)
        if not feature.dtype.startswith("cat"):
            expression = {"feature": feature, f"value{comparator}": value}
            feature_value = value_model.filter(**expression)
            new_expression["_feature_values__in"] = feature_value
        else:
            if isinstance(value, str):
                expression = {f"name{comparator}": value}
                label = ULabel.get(**expression)
                new_expression["ulabels"] = label
            else:
                raise NotImplementedError
    if cls == FeatureManager or cls == ParamManagerArtifact:
        return Artifact.filter(**new_expression)
    # might renable something similar in the future
    # elif cls == FeatureManagerCollection:
    #     return Collection.filter(**new_expression)
    elif cls == ParamManagerRun:
        return Run.filter(**new_expression)


@classmethod  # type: ignore
def filter(cls, **expression) -> QuerySet:
    """Query artifacts by features."""
    return filter_base(cls, **expression)


@classmethod  # type: ignore
def get(cls, **expression) -> Record:
    """Query a single artifact by feature."""
    return filter_base(cls, **expression).one()


@property  # type: ignore
def _feature_set_by_slot(self):
    """Feature sets by slot."""
    if self._feature_set_by_slot_ is None:
        self._feature_set_by_slot_ = get_feature_set_by_slot_(self._host)
    return self._feature_set_by_slot_


@property  # type: ignore
def _accessor_by_registry(self):
    """Accessor by ORM."""
    if self._accessor_by_registry_ is None:
        self._accessor_by_registry_ = get_accessor_by_registry_(self._host)
    return self._accessor_by_registry_


def _add_values(
    self,
    values: dict[str, str | int | float | bool],
    feature_param_field: FieldAttr,
    str_as_ulabel: bool = True,
) -> None:
    """Curate artifact with features & values.

    Args:
        values: A dictionary of keys (features) & values (labels, numbers, booleans).
        feature_param_field: The field of a reference registry to map keys of the
            dictionary.
    """
    # rename to distinguish from the values inside the dict
    features_values = values
    keys = features_values.keys()
    if isinstance(keys, DICT_KEYS_TYPE):
        keys = list(keys)  # type: ignore
    # deal with other cases later
    assert all(isinstance(key, str) for key in keys)  # noqa: S101
    registry = feature_param_field.field.model
    is_param = registry == Param
    model = Param if is_param else Feature
    value_model = ParamValue if is_param else FeatureValue
    model_name = "Param" if is_param else "Feature"
    if is_param:
        if self._host.__class__ == Artifact:
            if self._host.type != "model":
                raise ValidationError("Can only set params for model-like artifacts.")
    else:
        if self._host.__class__ == Artifact:
            if self._host.type != "dataset" and self._host.type is not None:
                raise ValidationError(
                    "Can only set features for dataset-like artifacts."
                )
    validated = registry.validate(keys, field=feature_param_field, mute=True)
    keys_array = np.array(keys)
    validated_keys = keys_array[validated]
    if validated.sum() != len(keys):
        not_validated_keys = keys_array[~validated]
        hint = "\n".join(
            [
                f"  ln.{model_name}(name='{key}', dtype='{infer_feature_type_convert_json(features_values[key], str_as_ulabel=str_as_ulabel)[0]}').save()"
                for key in not_validated_keys
            ]
        )
        msg = (
            f"These keys could not be validated: {not_validated_keys.tolist()}\n"
            f"Here is how to create a {model_name.lower()}:\n\n{hint}"
        )
        raise ValidationError(msg)
    registry.from_values(
        validated_keys,
        field=feature_param_field,
    )
    # figure out which of the values go where
    features_labels = defaultdict(list)
    _feature_values = []
    not_validated_values = []
    for key, value in features_values.items():
        feature = model.get(name=key)
        inferred_type, converted_value = infer_feature_type_convert_json(
            value,
            mute=True,
            str_as_ulabel=str_as_ulabel,
        )
        if feature.dtype == "number":
            if inferred_type not in {"int", "float"}:
                raise TypeError(
                    f"Value for feature '{key}' with type {feature.dtype} must be a number"
                )
        elif feature.dtype.startswith("cat"):
            if inferred_type != "?":
                if not (inferred_type.startswith("cat") or isinstance(value, Record)):
                    raise TypeError(
                        f"Value for feature '{key}' with type '{feature.dtype}' must be a string or record."
                    )
        elif not inferred_type == feature.dtype:
            raise ValidationError(
                f"Expected dtype for '{key}' is '{feature.dtype}', got '{inferred_type}'"
            )
        if not feature.dtype.startswith("cat"):
            # can remove the query once we have the unique constraint
            filter_kwargs = {model_name.lower(): feature, "value": converted_value}
            feature_value = value_model.filter(**filter_kwargs).one_or_none()
            if feature_value is None:
                feature_value = value_model(**filter_kwargs)
            _feature_values.append(feature_value)
        else:
            if isinstance(value, Record) or (
                isinstance(value, Iterable) and isinstance(next(iter(value)), Record)
            ):
                if isinstance(value, Record):
                    label_records = [value]
                else:
                    label_records = value  # type: ignore
                for record in label_records:
                    if record._state.adding:
                        raise ValidationError(
                            f"Please save {record} before annotation."
                        )
                    features_labels[record.__class__.__get_name_with_schema__()].append(
                        (feature, record)
                    )
            else:
                if isinstance(value, str):
                    values = [value]  # type: ignore
                else:
                    values = value  # type: ignore
                if "ULabel" not in feature.dtype:
                    feature.dtype += "[ULabel]"
                    feature.save()
                validated = ULabel.validate(values, field="name", mute=True)
                values_array = np.array(values)
                validated_values = values_array[validated]
                if validated.sum() != len(values):
                    not_validated_values += values_array[~validated].tolist()
                label_records = ULabel.from_values(validated_values, field="name")
                features_labels["ULabel"] += [
                    (feature, label_record) for label_record in label_records
                ]
    if not_validated_values:
        hint = (
            f"  ulabels = ln.ULabel.from_values({not_validated_values}, create=True)\n"
            f"  ln.save(ulabels)"
        )
        msg = (
            f"These values could not be validated: {not_validated_values}\n"
            f"Here is how to create ulabels for them:\n\n{hint}"
        )
        raise ValidationError(msg)
    # bulk add all links to ArtifactULabel
    if features_labels:
        if list(features_labels.keys()) != ["ULabel"]:
            related_names = dict_related_model_to_related_name(self._host.__class__)
        else:
            related_names = {"ULabel": "ulabels"}
        for class_name, registry_features_labels in features_labels.items():
            related_name = related_names[class_name]  # e.g., "ulabels"
            LinkORM = getattr(self._host, related_name).through
            field_name = f"{get_link_attr(LinkORM, self._host)}_id"  # e.g., ulabel_id
            links = [
                LinkORM(
                    **{
                        "artifact_id": self._host.id,
                        "feature_id": feature.id,
                        field_name: label.id,
                    }
                )
                for (feature, label) in registry_features_labels
            ]
            # a link might already exist
            try:
                save(links, ignore_conflicts=False)
            except Exception:
                save(links, ignore_conflicts=True)
                # now deal with links that were previously saved without a feature_id
                links_saved = LinkORM.filter(
                    **{
                        "artifact_id": self._host.id,
                        f"{field_name}__in": [
                            l.id for _, l in registry_features_labels
                        ],
                    }
                )
                for link in links_saved.all():
                    # TODO: also check for inconsistent features
                    if link.feature_id is None:
                        link.feature_id = [
                            f.id
                            for f, l in registry_features_labels
                            if l.id == getattr(link, field_name)
                        ][0]
                        link.save()
    if _feature_values:
        save(_feature_values)
        if is_param:
            LinkORM = self._host._param_values.through
            valuefield_id = "paramvalue_id"
        else:
            LinkORM = self._host._feature_values.through
            valuefield_id = "featurevalue_id"
        links = [
            LinkORM(
                **{
                    f"{self._host.__class__.__get_name_with_schema__().lower()}_id": self._host.id,
                    valuefield_id: feature_value.id,
                }
            )
            for feature_value in _feature_values
        ]
        # a link might already exist, to avoid raising a unique constraint
        # error, ignore_conflicts
        save(links, ignore_conflicts=True)


def add_values_features(
    self,
    values: dict[str, str | int | float | bool],
    feature_field: FieldAttr = Feature.name,
    str_as_ulabel: bool = True,
) -> None:
    """Curate artifact with features & values.

    Args:
        values: A dictionary of keys (features) & values (labels, numbers, booleans).
        feature_field: The field of a reference registry to map keys of the
            dictionary.
        str_as_ulabel: Whether to interpret string values as ulabels.
    """
    _add_values(self, values, feature_field, str_as_ulabel=str_as_ulabel)


def add_values_params(
    self,
    values: dict[str, str | int | float | bool],
) -> None:
    """Curate artifact with features & values.

    Args:
        values: A dictionary of keys (features) & values (labels, numbers, booleans).
    """
    _add_values(self, values, Param.name, str_as_ulabel=False)


def add_feature_set(self, feature_set: FeatureSet, slot: str) -> None:
    """Curate artifact with a feature set.

    Args:
        feature_set: `FeatureSet` A feature set record.
        slot: `str` The slot that marks where the feature set is stored in
            the artifact.
    """
    if self._host._state.adding:
        raise ValueError(
            "Please save the artifact or collection before adding a feature set!"
        )
    host_db = self._host._state.db
    feature_set.save(using=host_db)
    host_id_field = get_host_id_field(self._host)
    kwargs = {
        host_id_field: self._host.id,
        "featureset": feature_set,
        "slot": slot,
    }
    link_record = (
        self._host.feature_sets.through.objects.using(host_db)
        .filter(**kwargs)
        .one_or_none()
    )
    if link_record is None:
        self._host.feature_sets.through(**kwargs).save(using=host_db)
        if slot in self._feature_set_by_slot:
            logger.debug(f"replaced existing {slot} feature set")
        self._feature_set_by_slot_[slot] = feature_set  # type: ignore


def _add_set_from_df(
    self, field: FieldAttr = Feature.name, organism: str | None = None
):
    """Add feature set corresponding to column names of DataFrame."""
    if isinstance(self._host, Artifact):
        assert self._host._accessor == "DataFrame"  # noqa: S101
    else:
        # Collection
        assert self._host.artifact._accessor == "DataFrame"  # noqa: S101

    # parse and register features
    registry = field.field.model
    df = self._host.load()
    features = registry.from_values(df.columns, field=field, organism=organism)
    if len(features) == 0:
        logger.error(
            "no validated features found in DataFrame! please register features first!"
        )
        return

    # create and link feature sets
    feature_set = FeatureSet(features=features)
    feature_sets = {"columns": feature_set}
    self._host._feature_sets = feature_sets
    self._host.save()


def _add_set_from_anndata(
    self,
    var_field: FieldAttr,
    obs_field: FieldAttr | None = Feature.name,
    mute: bool = False,
    organism: str | Record | None = None,
):
    """Add features from AnnData."""
    if isinstance(self._host, Artifact):
        assert self._host._accessor == "AnnData"  # noqa: S101
    else:
        raise NotImplementedError()

    # parse and register features
    adata = self._host.load()
    feature_sets = parse_feature_sets_from_anndata(
        adata,
        var_field=var_field,
        obs_field=obs_field,
        mute=mute,
        organism=organism,
    )

    # link feature sets
    self._host._feature_sets = feature_sets
    self._host.save()


def _add_set_from_mudata(
    self,
    var_fields: dict[str, FieldAttr],
    obs_fields: dict[str, FieldAttr] = None,
    mute: bool = False,
    organism: str | Record | None = None,
):
    """Add features from MuData."""
    if obs_fields is None:
        obs_fields = {}
    if isinstance(self._host, Artifact):
        assert self._host._accessor == "MuData"  # noqa: S101
    else:
        raise NotImplementedError()

    # parse and register features
    mdata = self._host.load()
    feature_sets = {}
    obs_features = Feature.from_values(mdata.obs.columns)
    if len(obs_features) > 0:
        feature_sets["obs"] = FeatureSet(features=obs_features)
    for modality, field in var_fields.items():
        modality_fs = parse_feature_sets_from_anndata(
            mdata[modality],
            var_field=field,
            obs_field=obs_fields.get(modality, Feature.name),
            mute=mute,
            organism=organism,
        )
        for k, v in modality_fs.items():
            feature_sets[f"['{modality}'].{k}"] = v

    def unify_feature_sets_by_hash(feature_sets):
        unique_values = {}

        for key, value in feature_sets.items():
            value_hash = value.hash  # Assuming each value has a .hash attribute
            if value_hash in unique_values:
                feature_sets[key] = unique_values[value_hash]
            else:
                unique_values[value_hash] = value

        return feature_sets

    # link feature sets
    self._host._feature_sets = unify_feature_sets_by_hash(feature_sets)
    self._host.save()


def _add_from(self, data: Artifact | Collection, transfer_logs: dict = None):
    """Transfer features from a artifact or collection."""
    # This only covers feature sets
    if transfer_logs is None:
        transfer_logs = {"mapped": [], "transferred": []}
    using_key = settings._using_key
    for slot, feature_set in data.features._feature_set_by_slot.items():
        members = feature_set.members
        if len(members) == 0:
            continue
        registry = members[0].__class__
        # note here the features are transferred based on an unique field
        field = REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
        if hasattr(registry, "_ontology_id_field"):
            field = registry._ontology_id_field
        # this will be e.g. be a list of ontology_ids or uids
        member_uids = list(members.values_list(field, flat=True))
        # create records from ontology_id
        if hasattr(registry, "_ontology_id_field") and len(member_uids) > 0:
            # create from bionty
            members_records = registry.from_values(member_uids, field=field)
            save([r for r in members_records if r._state.adding])
        validated = registry.validate(member_uids, field=field, mute=True)
        new_members_uids = list(compress(member_uids, ~validated))
        new_members = members.filter(**{f"{field}__in": new_members_uids}).all()
        n_new_members = len(new_members)
        if n_new_members > 0:
            # transfer foreign keys needs to be run before transfer to default db
            transfer_fk_to_default_db_bulk(
                new_members, using_key, transfer_logs=transfer_logs
            )
            for feature in new_members:
                # not calling save=True here as in labels, because want to
                # bulk save below
                # transfer_fk is set to False because they are already transferred
                # in the previous step transfer_fk_to_default_db_bulk
                transfer_to_default_db(
                    feature, using_key, transfer_fk=False, transfer_logs=transfer_logs
                )
            logger.info(f"saving {n_new_members} new {registry.__name__} records")
            save(new_members)

        # create a new feature set from feature values using the same uid
        feature_set_self = FeatureSet.from_values(
            member_uids, field=getattr(registry, field)
        )
        if feature_set_self is None:
            if hasattr(registry, "organism_id"):
                logger.warning(
                    f"FeatureSet is not transferred, check if organism is set correctly: {feature_set}"
                )
            continue
        # make sure the uid matches if featureset is composed of same features
        if feature_set_self.hash == feature_set.hash:
            feature_set_self.uid = feature_set.uid
        logger.info(f"saving {slot} featureset: {feature_set_self}")
        self._host.features.add_feature_set(feature_set_self, slot)


FeatureManager.__init__ = __init__
ParamManager.__init__ = __init__
FeatureManager.__repr__ = __repr__
ParamManager.__repr__ = __repr__
FeatureManager.__getitem__ = __getitem__
FeatureManager.get_values = get_values
FeatureManager._feature_set_by_slot = _feature_set_by_slot
FeatureManager._accessor_by_registry = _accessor_by_registry
FeatureManager.add_values = add_values_features
FeatureManager.add_feature_set = add_feature_set
FeatureManager._add_set_from_df = _add_set_from_df
FeatureManager._add_set_from_anndata = _add_set_from_anndata
FeatureManager._add_set_from_mudata = _add_set_from_mudata
FeatureManager._add_from = _add_from
FeatureManager.filter = filter
FeatureManager.get = get
ParamManager.add_values = add_values_params
ParamManager.get_values = get_values
