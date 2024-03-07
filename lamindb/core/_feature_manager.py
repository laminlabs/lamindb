from itertools import compress
from typing import Dict, Iterable, Optional, Union

import anndata as ad
from anndata import AnnData
from lamin_utils import colors, logger
from lamindb_setup.core.upath import create_path
from lnschema_core.models import Artifact, Collection, Data, Feature, Registry
from lnschema_core.types import AnnDataLike, FieldAttr

from lamindb._feature import convert_numpy_dtype_to_lamin_feature_type
from lamindb._feature_set import FeatureSet
from lamindb._query_set import QuerySet
from lamindb._registry import (
    REGISTRY_UNIQUE_FIELD,
    get_default_str_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)
from lamindb._save import save
from lamindb.core.storage import LocalPathClasses

from ._settings import settings


def get_host_id_field(host: Union[Artifact, Collection]) -> str:
    if isinstance(host, Artifact):
        host_id_field = "artifact_id"
    else:
        host_id_field = "collection_id"
    return host_id_field


def get_accessor_by_orm(host: Union[Artifact, Collection]) -> Dict:
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
    host_db = host._state.db
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id}
    # otherwise, we need a query
    feature_set_links = host.feature_sets.through.objects.using(host_db).filter(
        **kwargs
    )
    return {
        feature_set_link.slot: FeatureSet.objects.using(host_db).get(
            id=feature_set_link.feature_set_id
        )
        for feature_set_link in feature_set_links
    }


def get_label_links(
    host: Union[Artifact, Collection], registry: str, feature: Feature
) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id, "feature_id": feature.id}
    link_records = (
        getattr(host, host.features._accessor_by_orm[registry])
        .through.objects.using(host._state.db)
        .filter(**kwargs)
    )
    return link_records


def get_feature_set_links(host: Union[Artifact, Collection]) -> QuerySet:
    host_id_field = get_host_id_field(host)
    kwargs = {host_id_field: host.id}
    feature_set_links = host.feature_sets.through.objects.filter(**kwargs)
    return feature_set_links


def print_features(self: Data) -> str:
    from lamindb._from_values import _print_values

    msg = ""
    features_lookup = Feature.objects.using(self._state.db).lookup().dict()
    for slot, feature_set in self.features._feature_set_by_slot.items():
        if feature_set.registry != "core.Feature":
            features = feature_set.members
            name_field = get_default_str_field(features[0])
            feature_names = [getattr(feature, name_field) for feature in features]
            msg += f"  {colors.bold(slot)}: {feature_set}\n"
            print_values = _print_values(feature_names, n=20)
            msg += f"    {print_values}\n"
        else:
            df_slot = feature_set.features.df()
            msg += f"  {colors.bold(slot)}: {feature_set}\n"
            for _, row in df_slot.iterrows():
                if row["type"] == "category" and row["registries"] is not None:
                    labels = self.labels.get(
                        features_lookup.get(row["name"]), mute=True
                    )
                    indent = ""
                    if isinstance(labels, dict):
                        msg += f"    🔗 {row['name']} ({row.registries})\n"
                        indent = "    "
                    else:
                        labels = {row["registries"]: labels}
                    for registry, labels in labels.items():  # noqa: B020
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


def parse_feature_sets_from_anndata(
    adata: AnnDataLike,
    var_field: FieldAttr,
    obs_field: FieldAttr = Feature.name,
    **kwargs,
) -> Dict:
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
        type = convert_numpy_dtype_to_lamin_feature_type(adata.X.dtype)
    feature_sets = {}
    logger.info("parsing feature names of X stored in slot 'var'")
    logger.indent = "   "
    feature_set_var = FeatureSet.from_values(
        data_parse.var.index,
        var_field,
        type=type,
        **kwargs,
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
            **kwargs,
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

    def __init__(self, host: Union[Artifact, Collection]):
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
        if hasattr(feature_set, "_features"):
            # feature set is not yet saved
            # need to think about turning this into a queryset
            return feature_set._features
        else:
            return getattr(feature_set, self._accessor_by_orm[orm_name]).all()

    def add(self, features: Iterable[Registry], slot: Optional[str] = None):
        """Add features stratified by slot."""
        if (hasattr(self._host, "accessor") and self._host.accessor == "DataFrame") or (
            hasattr(self._host, "artifact")
            and self._host.artifact.accessor == "DataFrame"
        ):
            slot = "columns" if slot is None else slot
        self._add_feature_set(feature_set=FeatureSet(features=features), slot=slot)

    def add_from_df(self):
        """Add features from DataFrame."""
        if isinstance(self._host, Artifact):
            assert self._host.accessor == "DataFrame"
        else:
            # Collection
            assert self._host.artifact.accessor == "DataFrame"

        # parse and register features
        df = self._host.load()
        features = Feature.from_values(df.columns)
        if len(features) == 0:
            logger.error(
                "no validated features found in DataFrame! please register features first:\n   → features = Feature.from_df(df)\n   → ln.save(features)"
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
        obs_field: Optional[FieldAttr] = Feature.name,
        **kwargs,
    ):
        """Add features from AnnData."""
        if isinstance(self._host, Artifact):
            assert self._host.accessor == "AnnData"
        else:
            # Collection
            assert self._host.artifact.accessor == "AnnData"

        # parse and register features
        adata = self._host.load()
        feature_sets = parse_feature_sets_from_anndata(
            adata, var_field=var_field, obs_field=obs_field, **kwargs
        )

        # link feature sets
        self._host._feature_sets = feature_sets
        self._host.save()

    def _add_feature_set(self, feature_set: FeatureSet, slot: str):
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
            self._feature_set_by_slot[slot] = feature_set

    def _add_from(self, data: Data, parents: bool = True):
        """Transfer features from a artifact or collection."""
        using_key = settings._using_key
        for slot, feature_set in data.features._feature_set_by_slot.items():
            members = feature_set.members
            if members.count() == 0:
                continue
            registry = members[0].__class__
            # note here the features are transferred based on an unique field
            field = REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
            if hasattr(registry, "ontology_id") and parents:
                field = "ontology_id"
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
            # TODO: make sure the uid matches if featureset is composed of same features
            # feature_set_self.uid = feature_set.uid
            logger.info(f"saving {slot} featureset: {feature_set_self}")
            self._host.features._add_feature_set(feature_set_self, slot)
