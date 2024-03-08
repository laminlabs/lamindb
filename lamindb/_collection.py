from collections import defaultdict
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Iterable,
    List,
    Literal,
    Optional,
    Tuple,
    Union,
)

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
from anndata import AnnData
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.hashing import hash_set
from lnschema_core.models import Collection, CollectionArtifact, FeatureSet
from lnschema_core.types import DataLike, VisibilityChoice

from lamindb._utils import attach_func_to_class_method
from lamindb.core._data import _track_run_input
from lamindb.core._mapped_collection import MappedCollection
from lamindb.core.versioning import get_uid_from_old_version, init_uid

from . import Artifact, Run
from ._artifact import data_is_anndata
from ._query_set import QuerySet
from ._registry import init_self_from_db
from .core._data import (
    add_transform_to_kwargs,
    get_run,
    save_feature_set_links,
    save_feature_sets,
)

if TYPE_CHECKING:
    from lamindb.core.storage._backed_access import AnnDataAccessor, BackedAccessor


def _check_accessor_collection(data: Any, accessor: Optional[str] = None):
    if accessor is None and isinstance(data, (AnnData, pd.DataFrame)):
        if isinstance(data, pd.DataFrame):
            logger.warning("data is a DataFrame, please use .from_df()")
            accessor = "DataFrame"
        elif data_is_anndata(data):
            logger.warning("data is an AnnData, please use .from_anndata()")
            accessor = "AnnData"
    return accessor


def __init__(
    collection: Collection,
    *args,
    **kwargs,
):
    if len(args) == len(collection._meta.concrete_fields):
        super(Collection, collection).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: data")
    data: Union[pd.DataFrame, ad.AnnData, Artifact, Iterable[Artifact]] = (
        kwargs.pop("data") if len(args) == 0 else args[0]
    )
    meta: Optional[str] = kwargs.pop("meta") if "meta" in kwargs else None
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    description: Optional[str] = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    reference: Optional[str] = (
        kwargs.pop("reference") if "reference" in kwargs else None
    )
    reference_type: Optional[str] = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    run: Optional[Run] = kwargs.pop("run") if "run" in kwargs else None
    is_new_version_of: Optional[Collection] = (
        kwargs.pop("is_new_version_of") if "is_new_version_of" in kwargs else None
    )
    version: Optional[str] = kwargs.pop("version") if "version" in kwargs else None
    visibility: Optional[int] = (
        kwargs.pop("visibility")
        if "visibility" in kwargs
        else VisibilityChoice.default.value
    )
    feature_sets: Dict[str, FeatureSet] = (
        kwargs.pop("feature_sets") if "feature_sets" in kwargs else {}
    )
    accessor = kwargs.pop("accessor") if "accessor" in kwargs else None
    if not isinstance(data, (Artifact, Iterable)):
        accessor = _check_accessor_collection(data=data, accessor=accessor)
    if not len(kwargs) == 0:
        raise ValueError(
            f"Only data, name, run, description, reference, reference_type, visibility can be passed, you passed: {kwargs}"
        )

    if is_new_version_of is None:
        provisional_uid = init_uid(version=version, n_full_id=20)
    else:
        if not isinstance(is_new_version_of, Collection):
            raise TypeError("is_new_version_of has to be of type ln.Collection")
        provisional_uid, version = get_uid_from_old_version(is_new_version_of, version)
        if name is None:
            name = is_new_version_of.name
    run = get_run(run)
    data_init_complete = False
    artifact = None
    artifacts = None
    # now handle potential metadata
    if meta is not None:
        if not isinstance(meta, (pd.DataFrame, ad.AnnData, Artifact)):
            raise ValueError(
                "meta has to be of type `(pd.DataFrame, ad.AnnData, Artifact)`"
            )
        data = meta
    # init artifact - is either data or metadata
    if isinstance(data, (pd.DataFrame, ad.AnnData, Artifact)):
        if isinstance(data, Artifact):
            artifact = data
            if artifact._state.adding:
                raise ValueError("Save artifact before creating collection!")
            if not feature_sets:
                feature_sets = artifact.features._feature_set_by_slot
            else:
                if len(artifact.features._feature_set_by_slot) > 0:
                    logger.info("overwriting feature sets linked to artifact")
        else:
            artifact_is_new_version_of = (
                is_new_version_of.artifact if is_new_version_of is not None else None
            )
            artifact = Artifact(
                data,
                run=run,
                description="tmp",
                version=version,
                is_new_version_of=artifact_is_new_version_of,
                accessor=accessor,
            )
            # do we really want to update the artifact here?
            if feature_sets:
                artifact._feature_sets = feature_sets
        hash = artifact.hash  # type: ignore
        provisional_uid = artifact.uid  # type: ignore
        if artifact.description is None or artifact.description == "tmp":
            artifact.description = f"See collection {provisional_uid}"  # type: ignore
        data_init_complete = True
    if not data_init_complete:
        if hasattr(data, "__getitem__"):
            assert isinstance(data[0], Artifact)  # type: ignore
            artifacts = data
            hash, feature_sets = from_artifacts(artifacts)  # type: ignore
            data_init_complete = True
        else:
            raise ValueError(
                "Only DataFrame, AnnData, Artifact or list of artifacts is allowed."
            )
    # we ignore collections in trash containing the same hash
    if hash is not None:
        existing_collection = Collection.filter(hash=hash).one_or_none()
    else:
        existing_collection = None
    if existing_collection is not None:
        logger.warning(
            f"returning existing collection with same hash: {existing_collection}"
        )
        init_self_from_db(collection, existing_collection)
        for slot, feature_set in collection.features._feature_set_by_slot.items():
            if slot in feature_sets:
                if not feature_sets[slot] == feature_set:
                    collection.feature_sets.remove(feature_set)
                    logger.warning(f"removing feature set: {feature_set}")
    else:
        kwargs = {}
        add_transform_to_kwargs(kwargs, run)
        super(Collection, collection).__init__(
            uid=provisional_uid,
            name=name,
            description=description,
            reference=reference,
            reference_type=reference_type,
            artifact=artifact,
            hash=hash,
            run=run,
            version=version,
            visibility=visibility,
            **kwargs,
        )
    collection._artifacts = artifacts
    collection._feature_sets = feature_sets
    # register provenance
    if is_new_version_of is not None:
        _track_run_input(is_new_version_of, run=run)
    if artifact is not None and artifact.run != run:
        _track_run_input(artifact, run=run)
    elif artifacts is not None:
        _track_run_input(artifacts, run=run)


@classmethod  # type: ignore
@doc_args(Collection.from_df.__doc__)
def from_df(
    cls,
    df: "pd.DataFrame",
    name: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
    reference: Optional[str] = None,
    reference_type: Optional[str] = None,
    version: Optional[str] = None,
    is_new_version_of: Optional["Artifact"] = None,
    **kwargs,
) -> "Collection":
    """{}."""
    if isinstance(df, Artifact):
        assert not df._state.adding
        assert df.accessor == "DataFrame"
    collection = Collection(
        data=df,
        name=name,
        run=run,
        description=description,
        reference=reference,
        reference_type=reference_type,
        version=version,
        is_new_version_of=is_new_version_of,
        accessor="DataFrame",
        **kwargs,
    )
    return collection


@classmethod  # type: ignore
@doc_args(Collection.from_anndata.__doc__)
def from_anndata(
    cls,
    adata: "AnnData",
    name: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
    reference: Optional[str] = None,
    reference_type: Optional[str] = None,
    version: Optional[str] = None,
    is_new_version_of: Optional["Artifact"] = None,
    **kwargs,
) -> "Collection":
    """{}."""
    if isinstance(adata, Artifact):
        assert not adata._state.adding
        assert adata.accessor == "AnnData"
    collection = Collection(
        data=adata,
        run=run,
        name=name,
        description=description,
        reference=reference,
        reference_type=reference_type,
        version=version,
        is_new_version_of=is_new_version_of,
        accessor="AnnData",
        **kwargs,
    )
    return collection


# internal function, not exposed to user
def from_artifacts(artifacts: Iterable[Artifact]) -> Tuple[str, Dict[str, str]]:
    # assert all artifacts are already saved
    logger.debug("check not saved")
    saved = not any(artifact._state.adding for artifact in artifacts)
    if not saved:
        raise ValueError("Not all artifacts are yet saved, please save them")
    # query all feature sets of artifacts
    logger.debug("artifact ids")
    artifact_ids = [artifact.id for artifact in artifacts]
    # query all feature sets at the same time rather
    # than making a single query per artifact
    logger.debug("feature_set_artifact_links")
    feature_set_artifact_links = Artifact.feature_sets.through.objects.filter(
        artifact_id__in=artifact_ids
    )
    feature_sets_by_slots = defaultdict(list)
    logger.debug("slots")
    for link in feature_set_artifact_links:
        feature_sets_by_slots[link.slot].append(link.feature_set_id)
    feature_sets_union = {}
    logger.debug("union")
    for slot, feature_set_ids_slot in feature_sets_by_slots.items():
        feature_set_1 = FeatureSet.filter(id=feature_set_ids_slot[0]).one()
        related_name = feature_set_1._get_related_name()
        features_registry = getattr(FeatureSet, related_name).field.model
        start_time = logger.debug("run filter")
        # this way of writing the __in statement turned out to be the fastest
        # evaluated on a link table with 16M entries connecting 500 feature sets with
        # 60k genes
        feature_ids = (
            features_registry.feature_sets.through.objects.filter(
                featureset_id__in=feature_set_ids_slot
            )
            .values(f"{features_registry.__name__.lower()}_id")
            .distinct()
        )
        start_time = logger.debug("done, start evaluate", time=start_time)
        features = features_registry.filter(id__in=feature_ids)
        feature_sets_union[slot] = FeatureSet(features, type=feature_set_1.type)
        start_time = logger.debug("done", time=start_time)
    # validate consistency of hashes
    # we do not allow duplicate hashes
    logger.debug("hashes")
    # artifact.hash is None for zarr
    # todo: more careful handling of such cases
    hashes = [artifact.hash for artifact in artifacts if artifact.hash is not None]
    if len(hashes) != len(set(hashes)):
        seen = set()
        non_unique = [x for x in hashes if x in seen or seen.add(x)]  # type: ignore
        raise ValueError(
            "Please pass artifacts with distinct hashes: these ones are non-unique"
            f" {non_unique}"
        )
    time = logger.debug("hash")
    hash = hash_set(set(hashes))
    logger.debug("done", time=time)
    return hash, feature_sets_union


# docstring handled through attach_func_to_class_method
def mapped(
    self,
    label_keys: Optional[Union[str, List[str]]] = None,
    join: Optional[Literal["inner", "outer"]] = "inner",
    encode_labels: Union[bool, List[str]] = True,
    unknown_label: Optional[Union[str, Dict[str, str]]] = None,
    cache_categories: bool = True,
    parallel: bool = False,
    dtype: Optional[str] = None,
    stream: bool = False,
    is_run_input: Optional[bool] = None,
) -> "MappedCollection":
    _track_run_input(self, is_run_input)
    path_list = []
    for artifact in self.artifacts.all():
        if artifact.suffix not in {".h5ad", ".zrad", ".zarr"}:
            logger.warning(f"Ignoring artifact with suffix {artifact.suffix}")
            continue
        elif not stream and artifact.suffix == ".h5ad":
            path_list.append(artifact.stage())
        else:
            path_list.append(artifact.path)
    return MappedCollection(
        path_list,
        label_keys,
        join,
        encode_labels,
        unknown_label,
        cache_categories,
        parallel,
        dtype,
    )


# docstring handled through attach_func_to_class_method
def backed(
    self, is_run_input: Optional[bool] = None
) -> Union["AnnDataAccessor", "BackedAccessor"]:
    _track_run_input(self, is_run_input)
    if self.artifact is None:
        raise RuntimeError(
            "Can only call backed() for collections with a single artifact"
        )
    return self.artifact.backed()


# docstring handled through attach_func_to_class_method
def load(
    self,
    join: Literal["inner", "outer"] = "outer",
    is_run_input: Optional[bool] = None,
    **kwargs,
) -> DataLike:
    # cannot call _track_run_input here, see comment further down
    if self.artifact is not None:
        _track_run_input(self, is_run_input)
        return self.artifact.load()
    else:
        all_artifacts = self.artifacts.all()
        suffixes = [artifact.suffix for artifact in all_artifacts]
        if len(set(suffixes)) != 1:
            raise RuntimeError(
                "Can only load collections where all artifacts have the same suffix"
            )
        # because we're tracking data flow on the collection-level, here, we don't
        # want to track it on the artifact-level
        objects = [artifact.load(is_run_input=False) for artifact in all_artifacts]
        artifact_uids = [artifact.uid for artifact in all_artifacts]
        if isinstance(objects[0], pd.DataFrame):
            concat_object = pd.concat(objects, join=join)
        elif isinstance(objects[0], ad.AnnData):
            concat_object = ad.concat(
                objects, join=join, label="artifact_uid", keys=artifact_uids
            )
        # only call it here because there might be errors during concat
        _track_run_input(self, is_run_input)
        return concat_object


# docstring handled through attach_func_to_class_method
def delete(self, permanent: Optional[bool] = None) -> None:
    # change visibility to trash
    trash_visibility = VisibilityChoice.trash.value
    if self.visibility > trash_visibility and permanent is not True:
        self.visibility = trash_visibility
        self.save()
        logger.warning(f"moved collection to trash (visibility = {trash_visibility})")
        return

    # permanent delete
    if permanent is None:
        response = input(
            "Collection record is already in trash! Are you sure to delete it from your"
            " database? (y/n) You can't undo this action."
        )
        delete_record = response == "y"
    else:
        delete_record = permanent

    if delete_record:
        super(Collection, self).delete()


# docstring handled through attach_func_to_class_method
def save(self, *args, **kwargs) -> None:
    if self.artifact is not None:
        self.artifact.save()
    # we don't need to save feature sets again
    save_feature_sets(self)
    super(Collection, self).save()
    # we don't allow updating the collection of artifacts
    # if users want to update the set of artifacts, they
    # have to create a new collection
    if hasattr(self, "_artifacts"):
        if self._artifacts is not None and len(self._artifacts) > 0:
            links = [
                CollectionArtifact(collection_id=self.id, artifact_id=artifact.id)
                for artifact in self._artifacts
            ]
            # the below seems to preserve the order of the list in the
            # auto-incrementing integer primary
            # merely using .unordered_artifacts.set(*...) doesn't achieve this
            # we need ignore_conflicts=True so that this won't error if links already exist
            CollectionArtifact.objects.bulk_create(links, ignore_conflicts=True)
    save_feature_set_links(self)


# docstring handled through attach_func_to_class_method
def restore(self) -> None:
    self.visibility = VisibilityChoice.default.value
    self.save()
    if self.artifact is not None:
        self.artifact.visibility = VisibilityChoice.default.value
        self.artifact.save()


@property  # type: ignore
@doc_args(Collection.artifacts.__doc__)
def artifacts(self) -> QuerySet:
    """{}."""
    _track_run_input(self)
    return self.unordered_artifacts.order_by("collectionartifact__id")


METHOD_NAMES = [
    "__init__",
    "from_anndata",
    "from_df",
    "mapped",
    "backed",
    "load",
    "delete",
    "save",
    "restore",
]

if ln_setup._TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(Collection, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Collection, globals())

# this seems a Django-generated function
delattr(Collection, "get_visibility_display")
Collection.artifacts = artifacts
