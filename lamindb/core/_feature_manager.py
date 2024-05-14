from __future__ import annotations

from itertools import compress
from typing import TYPE_CHECKING, Iterable, Optional

import anndata as ad
from anndata import AnnData
from lamin_utils import colors, logger
from lamindb_setup.core.upath import create_path
from lnschema_core.models import Artifact, Collection, Data, Feature, Registry

from lamindb._feature import convert_numpy_dtype_to_lamin_feature_type
from lamindb._feature_set import FeatureSet
from lamindb._registry import (
    REGISTRY_UNIQUE_FIELD,
    get_default_str_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)
from lamindb._save import save
from lamindb.core.storage import LocalPathClasses

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


def get_accessor_by_orm(host: Artifact | Collection) -> dict:
    dictionary = {
        field.related_model.__get_name_with_schema__(): field.name
        for field in host._meta.related_objects
    }
    dictionary["core.Feature"] = "features"
    dictionary["core.ULabel"] = "ulabels"
    return dictionary


def get_feature_set_by_slot(host) -> dict:
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
        .select_related("feature_set")
    )

    return {fsl.slot: fsl.feature_set for fsl in feature_set_links}


def get_label_links(
    host: Artifact | Collection, registry: str, feature: Feature
) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id, "feature_id": feature.id}
    link_records = (
        getattr(host, host.features.accessor_by_orm[registry])
        .through.objects.using(host._state.db)
        .filter(**kwargs)
    )
    return link_records


def get_feature_set_links(host: Artifact | Collection) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id}
    feature_set_links = host.feature_sets.through.objects.filter(**kwargs)
    return feature_set_links


def print_features(self: Data) -> str:
    from lamindb._from_values import _print_values

    from ._data import format_repr

    messages = []
    for slot, feature_set in get_feature_set_by_slot(self).items():
        if feature_set.registry != "core.Feature":
            features = feature_set.members
            # features.first() is a lot slower than features[0] here
            name_field = get_default_str_field(features[0])
            feature_names = list(features.values_list(name_field, flat=True)[:30])
            messages.append(
                f"  {colors.bold(slot)}: {format_repr(feature_set, exclude='hash')}\n"
            )
            print_values = _print_values(feature_names, n=20)
            messages.append(f"    {print_values}\n")
        else:
            features_lookup = {
                f.name: f for f in Feature.objects.using(self._state.db).filter().all()
            }
            messages.append(
                f"  {colors.bold(slot)}: {format_repr(feature_set, exclude='hash')}\n"
            )
            for name, row_type, registries in feature_set.features.values_list(
                "name", "type", "registries"
            ):
                if row_type == "category" and registries is not None:
                    labels = self.labels.get(features_lookup.get(name), mute=True)
                    indent = ""
                    if isinstance(labels, dict):
                        messages.append(f"    🔗 {name} ({registries})\n")
                        indent = "    "
                    else:
                        labels = {registries: labels}
                    for registry, registry_labels in labels.items():
                        field = get_default_str_field(registry_labels)
                        values_list = registry_labels.values_list(field, flat=True)
                        count_str = f"{feature_set.n}, {colors.italic(f'{registry}')}"
                        print_values = _print_values(values_list[:20], n=10)
                        msg_objects = (
                            f"{indent}    🔗 {name} ({count_str}):" f" {print_values}\n"
                        )
                        messages.append(msg_objects)
                else:
                    messages.append(f"    {name} ({row_type})\n")
    if messages:
        messages.insert(0, f"{colors.green('Features')}:\n")
    return "".join(messages)


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


class FeatureManager:
    """Feature manager (:attr:`~lamindb.core.Data.features`).

    See :class:`~lamindb.core.Data` for more information.
    """

    def __init__(self, host: Artifact | Collection):
        self._host = host
        self._feature_set_by_slot = None
        self._accessor_by_orm = None

    def __repr__(self) -> str:
        if len(self.feature_set_by_slot) > 0:
            return print_features(self._host)
        else:
            return "no linked features"

    def __getitem__(self, slot) -> QuerySet:
        if slot not in self.feature_set_by_slot:
            raise ValueError(
                f"No linked feature set for slot: {slot}\nDid you get validation"
                " warnings? Only features that match registered features get validated"
                " and linked."
            )
        feature_set = self.feature_set_by_slot[slot]
        orm_name = feature_set.registry
        if hasattr(feature_set, "_features"):
            # feature set is not yet saved
            # need to think about turning this into a queryset
            return feature_set._features
        else:
            return getattr(feature_set, self.accessor_by_orm[orm_name]).all()

    @property
    def feature_set_by_slot(self):
        """Feature sets by slot."""
        if self._feature_set_by_slot is None:
            self._feature_set_by_slot = get_feature_set_by_slot(self._host)
        return self._feature_set_by_slot

    @property
    def accessor_by_orm(self):
        """Accessor by ORM."""
        if self._accessor_by_orm is None:
            self._accessor_by_orm = get_accessor_by_orm(self._host)
        return self._accessor_by_orm

    def add(self, features: Iterable[Registry], slot: str | None = None):
        """Add features stratified by slot."""
        if (hasattr(self._host, "accessor") and self._host.accessor == "DataFrame") or (
            hasattr(self._host, "artifact")
            and self._host.artifact.accessor == "DataFrame"
        ):
            slot = "columns" if slot is None else slot
        self.add_feature_set(feature_set=FeatureSet(features=features), slot=slot)

    def add_from_df(self, field: FieldAttr = Feature.name, organism: str | None = None):
        """Add features from DataFrame."""
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

    def add_from_anndata(
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

    def add_from_mudata(
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

    def add_feature_set(self, feature_set: FeatureSet, slot: str):
        """Add new feature set to a slot.

        Args:
            feature_set: `FeatureSet` A feature set object.
            slot: `str` The access slot.
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
            "feature_set": feature_set,
            "slot": slot,
        }
        link_record = (
            self._host.feature_sets.through.objects.using(host_db)
            .filter(**kwargs)
            .one_or_none()
        )
        if link_record is None:
            self._host.feature_sets.through(**kwargs).save(using=host_db)
            if slot in self.feature_set_by_slot:
                logger.warning(f"replaced existing {slot} featureset")
            # this _feature_set_by_slot here is private
            self._feature_set_by_slot[slot] = feature_set  # type: ignore

    def _add_from(self, data: Data, parents: bool = True):
        """Transfer features from a artifact or collection."""
        using_key = settings._using_key
        for slot, feature_set in data.features.feature_set_by_slot.items():
            members = feature_set.members
            if members.count() == 0:
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
            if field == "ontology_id" and len(member_uids) > 0:
                # create from bionty
                records = registry.from_values(member_uids, field=field)
                if len(records) > 0:
                    save(records, parents=parents)
            validated = registry.validate(member_uids, field=field, mute=True)
            new_members_uids = list(compress(member_uids, ~validated))
            new_members = members.filter(**{f"{field}__in": new_members_uids}).all()
            if new_members.count() > 0:
                mute = True if new_members.count() > 10 else False
                # transfer foreign keys needs to be run before transfer to default db
                transfer_fk_to_default_db_bulk(new_members, using_key)
                for feature in new_members:
                    # not calling save=True here as in labels, because want to
                    # bulk save below
                    # transfer_fk is set to False because they are already transferred
                    # in the previous step transfer_fk_to_default_db_bulk
                    transfer_to_default_db(
                        feature, using_key, mute=mute, transfer_fk=False
                    )
                logger.info(
                    f"saving {new_members.count()} new {registry.__name__} records"
                )
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
