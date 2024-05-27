from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from itertools import compress
from typing import TYPE_CHECKING, Any

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData
from lamin_utils import colors, logger
from lamindb_setup.core.upath import create_path
from lnschema_core.models import (
    Artifact,
    Collection,
    Data,
    Feature,
    FeatureManager,
    FeatureValue,
    LinkORM,
    Registry,
    ULabel,
)

from lamindb._feature import FEATURE_TYPES, convert_numpy_dtype_to_lamin_feature_type
from lamindb._feature_set import DICT_KEYS_TYPE, FeatureSet
from lamindb._registry import (
    REGISTRY_UNIQUE_FIELD,
    get_default_str_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)
from lamindb._save import save
from lamindb.core.exceptions import ValidationError
from lamindb.core.storage import LocalPathClasses

from ._label_manager import get_labels_as_dict
from ._settings import settings

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
    feature_set_links = (
        host.feature_sets.through.objects.using(host_db)
        .filter(**kwargs)
        .select_related("featureset")
    )
    return {fsl.slot: fsl.featureset for fsl in feature_set_links}


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
    feature_set_links = host.feature_sets.through.objects.filter(**kwargs)
    return feature_set_links


def get_link_attr(link: LinkORM, data: Data) -> str:
    link_model_name = link.__class__.__name__
    link_attr = link_model_name.replace(data.__class__.__name__, "")
    if link_attr == "ExperimentalFactor":
        link_attr = "experimental_factor"
    else:
        link_attr = link_attr.lower()
    return link_attr


def print_features(self: Data, print_types: bool = False) -> str:
    from lamindb._from_values import _print_values

    msg = ""
    # feature values
    labels_msg = ""
    labels_by_feature = defaultdict(list)
    for _, (_, links) in get_labels_as_dict(self, links=True).items():
        for link in links:
            if link.feature_id is not None:
                link_attr = get_link_attr(link, self)
                labels_by_feature[link.feature_id].append(getattr(link, link_attr).name)
    for feature_id, labels_list in labels_by_feature.items():
        feature = Feature.objects.using(self._state.db).get(id=feature_id)
        print_values = _print_values(labels_list, n=10)
        type_str = f": {feature.dtype}" if print_types else ""
        labels_msg += f"    '{feature.name}'{type_str} = {print_values}\n"
    if labels_msg:
        msg += f"  {colors.italic('Features')}\n"
        msg += labels_msg
    # feature sets
    feature_set_msg = ""
    for slot, feature_set in get_feature_set_by_slot_(self).items():
        features = feature_set.members
        # features.first() is a lot slower than features[0] here
        name_field = get_default_str_field(features[0])
        feature_names = list(features.values_list(name_field, flat=True)[:20])
        type_str = f": {feature_set.registry}" if print_types else ""
        feature_set_msg += f"    '{slot}'{type_str} = {_print_values(feature_names)}\n"
    if feature_set_msg:
        msg += f"  {colors.italic('Feature sets')}\n"
        msg += feature_set_msg
    return msg


def parse_feature_sets_from_anndata(
    adata: AnnData,
    var_field: FieldAttr | None = None,
    obs_field: FieldAttr = Feature.name,
    mute: bool = False,
    organism: str | Registry | None = None,
) -> dict:
    data_parse = adata
    if not isinstance(adata, AnnData):  # is a path
        filepath = create_path(adata)  # returns Path for local
        if not isinstance(filepath, LocalPathClasses):
            from lamindb.core.storage._backed_access import backed_access

            using_key = settings._using_key
            data_parse = backed_access(filepath, using_key)
        else:
            data_parse = ad.read(filepath, backed="r")
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


def infer_feature_type_convert_json(value: Any, mute: bool = False) -> tuple[str, Any]:
    if isinstance(value, bool):
        return FEATURE_TYPES["bool"], value
    elif isinstance(value, int):
        return FEATURE_TYPES["int"], value
    elif isinstance(value, float):
        return FEATURE_TYPES["float"], value
    elif isinstance(value, str):
        return FEATURE_TYPES["str"] + "[ULabel]", value
    elif isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        if isinstance(value, (pd.Series, np.ndarray)):
            return convert_numpy_dtype_to_lamin_feature_type(value.dtype), list(value)
        if len(value) > 0:  # type: ignore
            first_element_type = type(next(iter(value)))
            if all(isinstance(elem, first_element_type) for elem in value):
                if first_element_type == bool:
                    return FEATURE_TYPES["bool"], value
                elif first_element_type == int:
                    return FEATURE_TYPES["int"], value
                elif first_element_type == float:
                    return FEATURE_TYPES["float"], value
                elif first_element_type == str:
                    return FEATURE_TYPES["str"] + "[ULabel]", value
    if not mute:
        logger.warning(f"cannot infer feature type of: {value}, returning '?")
    return ("?", value)


def __init__(self, host: Artifact | Collection):
    self._host = host
    self._feature_set_by_slot_ = None
    self._accessor_by_registry_ = None


def __repr__(self) -> str:
    return print_features(self._host)


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


@classmethod  # type:ignore
def filter(cls) -> QuerySet:
    """Filter features."""
    print("cls.__name__", cls.__name__)


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


def add_values(
    self,
    values: dict[str, str | int | float | bool],
    feature_field: FieldAttr = Feature.name,
) -> None:
    """Annotate artifact with features & values.

    Args:
        values: A dictionary of keys (features) & values (labels, numbers, booleans).
        feature_field: The field of a reference registry to map keys of the
            dictionary.
    """
    # rename to distinguish from the values inside the dict
    features_values = values
    keys = features_values.keys()
    if isinstance(keys, DICT_KEYS_TYPE):
        keys = list(keys)  # type: ignore
    # deal with other cases later
    assert all(isinstance(key, str) for key in keys)
    registry = feature_field.field.model
    validated = registry.validate(keys, field=feature_field, mute=True)
    keys_array = np.array(keys)
    validated_keys = keys_array[validated]
    if validated.sum() != len(keys):
        not_validated_keys = keys_array[~validated]
        hint = "\n".join(
            [
                f"  ln.Feature(name='{key}', dtype='{infer_feature_type_convert_json(features_values[key])[0]}').save()"
                for key in not_validated_keys
            ]
        )
        msg = (
            f"These keys could not be validated: {not_validated_keys.tolist()}\n"
            f"If there are no typos, create features for them:\n\n{hint}"
        )
        raise ValidationError(msg)
    registry.from_values(
        validated_keys,
        field=feature_field,
    )
    # figure out which of the values go where
    features_labels = []
    feature_values = []
    not_validated_values = []
    for key, value in features_values.items():
        feature = Feature.filter(name=key).one()
        inferred_type, converted_value = infer_feature_type_convert_json(
            value, mute=True
        )
        if feature.dtype == "number":
            if inferred_type not in {"int", "float"}:
                raise TypeError(
                    f"Value for feature '{key}' with type {feature.dtype} must be a number"
                )
        elif feature.dtype == "cat":
            if not (inferred_type.startswith("cat") or isinstance(value, Registry)):
                raise TypeError(
                    f"Value for feature '{key}' with type '{feature.dtype}' must be a string or record."
                )
        elif feature.dtype == "bool":
            assert isinstance(value, bool)
        if not feature.dtype.startswith("cat"):
            feature_values.append(FeatureValue(feature=feature, value=converted_value))
        else:
            if isinstance(value, Registry):
                assert not value._state.adding
                label_record = value
                assert isinstance(label_record, ULabel)
                features_labels.append((feature, label_record))
            else:
                if isinstance(value, str):
                    values = [value]  # type: ignore
                else:
                    values = value  # type: ignore
                validated = ULabel.validate(values, field="name", mute=True)
                values_array = np.array(values)
                validated_values = values_array[validated]
                if validated.sum() != len(values):
                    not_validated_values += values_array[~validated].tolist()
                label_records = ULabel.from_values(validated_values, field="name")
                features_labels += [
                    (feature, label_record) for label_record in label_records
                ]
    if not_validated_values:
        hint = (
            f"  ulabels = ln.ULabel.from_values({not_validated_values}, create=True)\n"
            f"  ln.save(ulabels)"
        )
        msg = (
            f"These values could not be validated: {not_validated_values}\n"
            f"If there are no typos, create ulabels for them:\n\n{hint}"
        )
        raise ValidationError(msg)
    # bulk add all links to ArtifactULabel
    if features_labels:
        LinkORM = self._host.ulabels.through
        links = [
            LinkORM(
                artifact_id=self._host.id, feature_id=feature.id, ulabel_id=label.id
            )
            for (feature, label) in features_labels
        ]
        # a link might already exist, to avoid raising a unique constraint
        # error, ignore_conflicts
        save(links, ignore_conflicts=True)
    if feature_values:
        save(feature_values)
        LinkORM = self._host.feature_values.through
        links = [
            LinkORM(artifact_id=self._host.id, featurevalue_id=feature_value.id)
            for feature_value in feature_values
        ]
        save(links)


def add_feature_set(self, feature_set: FeatureSet, slot: str) -> None:
    """Annotate artifact with a feature set.

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
        assert self._host.accessor == "DataFrame"
    else:
        # Collection
        assert self._host.artifact.accessor == "DataFrame"

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
    organism: str | Registry | None = None,
):
    """Add features from AnnData."""
    if isinstance(self._host, Artifact):
        assert self._host.accessor == "AnnData"
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
    organism: str | Registry | None = None,
):
    """Add features from MuData."""
    if obs_fields is None:
        obs_fields = {}
    if isinstance(self._host, Artifact):
        assert self._host.accessor == "MuData"
    else:
        raise NotImplementedError()

    # parse and register features
    mdata = self._host.load()
    feature_sets = {}
    obs_features = features = Feature.from_values(mdata.obs.columns)
    if len(obs_features) > 0:
        feature_sets["obs"] = FeatureSet(features=features)
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

    # link feature sets
    self._host._feature_sets = feature_sets
    self._host.save()


def _add_from(self, data: Data, parents: bool = True):
    """Transfer features from a artifact or collection."""
    # This only covers feature sets, though.
    using_key = settings._using_key
    for slot, feature_set in data.features._feature_set_by_slot.items():
        members = feature_set.members
        if len(members) == 0:
            continue
        registry = members[0].__class__
        # note here the features are transferred based on an unique field
        field = REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
        # TODO: get a default ID field for the registry
        if hasattr(registry, "ontology_id") and parents:
            field = "ontology_id"
        elif hasattr(registry, "ensembl_gene_id"):
            field = "ensembl_gene_id"
        elif hasattr(registry, "uniprotkb_id"):
            field = "uniprotkb_id"

        if registry.__get_name_with_schema__() == "bionty.Organism":
            parents = False
        # this will be e.g. be a list of ontology_ids or uids
        member_uids = list(members.values_list(field, flat=True))
        # create records from ontology_id in order to populate parents
        if field == "ontology_id" and len(member_uids) > 0 and parents:
            # create from bionty
            records = registry.from_values(member_uids, field=field)
            if len(records) > 0:
                save(records, parents=parents)
        validated = registry.validate(member_uids, field=field, mute=True)
        new_members_uids = list(compress(member_uids, ~validated))
        new_members = members.filter(**{f"{field}__in": new_members_uids}).all()
        n_new_members = len(new_members)
        if n_new_members > 0:
            mute = True if n_new_members > 10 else False
            # transfer foreign keys needs to be run before transfer to default db
            transfer_fk_to_default_db_bulk(new_members, using_key)
            for feature in new_members:
                # not calling save=True here as in labels, because want to
                # bulk save below
                # transfer_fk is set to False because they are already transferred
                # in the previous step transfer_fk_to_default_db_bulk
                transfer_to_default_db(feature, using_key, mute=mute, transfer_fk=False)
            logger.info(f"saving {n_new_members} new {registry.__name__} records")
            save(new_members, parents=parents)

        # create a new feature set from feature values using the same uid
        feature_set_self = FeatureSet.from_values(
            member_uids, field=getattr(registry, field)
        )
        if feature_set_self is None:
            if hasattr(registry, "organism"):
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
FeatureManager.__repr__ = __repr__
FeatureManager.__getitem__ = __getitem__
FeatureManager._feature_set_by_slot = _feature_set_by_slot
FeatureManager._accessor_by_registry = _accessor_by_registry
FeatureManager.add_values = add_values
FeatureManager.add_feature_set = add_feature_set
FeatureManager._add_set_from_df = _add_set_from_df
FeatureManager._add_set_from_anndata = _add_set_from_anndata
FeatureManager._add_set_from_mudata = _add_set_from_mudata
FeatureManager._add_from = _add_from
FeatureManager.filter = filter
