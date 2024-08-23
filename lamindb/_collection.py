from __future__ import annotations

from collections import defaultdict
from typing import (
    TYPE_CHECKING,
    Any,
    Iterable,
    Literal,
)

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.hashing import hash_set
from lnschema_core.models import (
    Collection,
    CollectionArtifact,
    FeatureManager,
    FeatureSet,
)
from lnschema_core.types import VisibilityChoice

from lamindb._artifact import update_attributes
from lamindb._utils import attach_func_to_class_method
from lamindb.core._data import _track_run_input
from lamindb.core._mapped_collection import MappedCollection
from lamindb.core.versioning import get_uid_from_old_version, init_uid, process_revises

from . import Artifact, Run
from ._record import init_self_from_db
from .core._data import (
    add_transform_to_kwargs,
    get_run,
    save_feature_set_links,
    save_feature_sets,
)
from .core._settings import settings

if TYPE_CHECKING:
    from lamindb.core.storage import UPath

    from ._query_set import QuerySet


def __init__(
    collection: Collection,
    *args,
    **kwargs,
):
    collection.features = FeatureManager(collection)
    if len(args) == len(collection._meta.concrete_fields):
        super(Collection, collection).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: artifacts")
    artifacts: Artifact | Iterable[Artifact] = (
        kwargs.pop("artifacts") if len(args) == 0 else args[0]
    )
    meta_artifact: Artifact | None = (
        kwargs.pop("meta_artifact") if "meta_artifact" in kwargs else None
    )
    name: str | None = kwargs.pop("name") if "name" in kwargs else None
    description: str | None = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    reference: str | None = kwargs.pop("reference") if "reference" in kwargs else None
    reference_type: str | None = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    run: Run | None = kwargs.pop("run") if "run" in kwargs else None
    revises: Collection | None = kwargs.pop("revises") if "revises" in kwargs else None
    version: str | None = kwargs.pop("version") if "version" in kwargs else None
    visibility: int | None = (
        kwargs.pop("visibility")
        if "visibility" in kwargs
        else VisibilityChoice.default.value
    )
    feature_sets: dict[str, FeatureSet] = (
        kwargs.pop("feature_sets") if "feature_sets" in kwargs else {}
    )
    if not len(kwargs) == 0:
        raise ValueError(
            f"Only artifacts, name, run, description, reference, reference_type, visibility can be passed, you passed: {kwargs}"
        )
    provisional_uid, version, name = process_revises(revises, version, name, Collection)
    run = get_run(run)
    if isinstance(artifacts, Artifact):
        artifacts = [artifacts]
    else:
        if not hasattr(artifacts, "__getitem__"):
            raise ValueError("Artifact or List[Artifact] is allowed.")
        assert isinstance(artifacts[0], Artifact)  # type: ignore  # noqa: S101
    hash, feature_sets = from_artifacts(artifacts)  # type: ignore
    if meta_artifact is not None:
        if not isinstance(meta_artifact, Artifact):
            raise ValueError("meta_artifact has to be an Artifact")
        if isinstance(meta_artifact, Artifact):
            if meta_artifact._state.adding:
                raise ValueError(
                    "Save meta_artifact artifact before creating collection!"
                )
            if not feature_sets:
                feature_sets = meta_artifact.features._feature_set_by_slot
            else:
                if len(meta_artifact.features._feature_set_by_slot) > 0:
                    logger.info("overwriting feature sets linked to artifact")
    # we ignore collections in trash containing the same hash
    if hash is not None:
        existing_collection = Collection.filter(hash=hash).one_or_none()
    else:
        existing_collection = None
    if existing_collection is not None:
        logger.warning(
            f"returning existing collection with same hash: {existing_collection}"
        )
        # update the run of the existing artifact
        if run is not None:
            # save the information that this artifact was previously
            # produced by another run
            if existing_collection.run is not None:
                existing_collection.run.output_collections_with_later_updates.add(
                    existing_collection
                )
            # update the run of the artifact with the latest run
            existing_collection.run = run
            existing_collection.transform = run.transform
        init_self_from_db(collection, existing_collection)
        update_attributes(collection, {"description": description, "name": name})
        for slot, feature_set in collection.features._feature_set_by_slot.items():
            if slot in feature_sets:
                if not feature_sets[slot] == feature_set:
                    collection.feature_sets.remove(feature_set)
                    logger.warning(f"removing feature set: {feature_set}")
    else:
        kwargs = {}
        add_transform_to_kwargs(kwargs, run)
        search_names_setting = settings.creation.search_names
        if revises is not None and name == revises.name:
            settings.creation.search_names = False
        super(Collection, collection).__init__(
            uid=provisional_uid,
            name=name,
            description=description,
            reference=reference,
            reference_type=reference_type,
            meta_artifact=meta_artifact,
            hash=hash,
            run=run,
            version=version,
            visibility=visibility,
            revises=revises,
            **kwargs,
        )
        settings.creation.search_names = search_names_setting
    collection._artifacts = artifacts
    collection._feature_sets = feature_sets
    # register provenance
    if revises is not None:
        _track_run_input(revises, run=run)
    _track_run_input(artifacts, run=run)


# internal function, not exposed to user
def from_artifacts(artifacts: Iterable[Artifact]) -> tuple[str, dict[str, str]]:
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
    logger.debug("links_feature_set_artifact")
    links_feature_set_artifact = Artifact.feature_sets.through.objects.filter(
        artifact_id__in=artifact_ids
    )
    feature_sets_by_slots = defaultdict(list)
    logger.debug("slots")
    for link in links_feature_set_artifact:
        feature_sets_by_slots[link.slot].append(link.featureset_id)
    feature_sets_union = {}
    logger.debug("union")
    for slot, feature_set_ids_slot in feature_sets_by_slots.items():
        feature_set_1 = FeatureSet.get(id=feature_set_ids_slot[0])
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
        feature_sets_union[slot] = FeatureSet(features, dtype=feature_set_1.dtype)
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
    layers_keys: str | list[str] | None = None,
    obs_keys: str | list[str] | None = None,
    obsm_keys: str | list[str] | None = None,
    join: Literal["inner", "outer"] | None = "inner",
    encode_labels: bool | list[str] = True,
    unknown_label: str | dict[str, str] | None = None,
    cache_categories: bool = True,
    parallel: bool = False,
    dtype: str | None = None,
    stream: bool = False,
    is_run_input: bool | None = None,
) -> MappedCollection:
    path_list = []
    for artifact in self.ordered_artifacts.all():
        if artifact.suffix not in {".h5ad", ".zarr"}:
            logger.warning(f"Ignoring artifact with suffix {artifact.suffix}")
            continue
        elif not stream:
            path_list.append(artifact.cache())
        else:
            path_list.append(artifact.path)
    ds = MappedCollection(
        path_list,
        layers_keys,
        obs_keys,
        obsm_keys,
        join,
        encode_labels,
        unknown_label,
        cache_categories,
        parallel,
        dtype,
    )
    # track only if successful
    _track_run_input(self, is_run_input)
    return ds


# docstring handled through attach_func_to_class_method
def cache(self, is_run_input: bool | None = None) -> list[UPath]:
    path_list = []
    for artifact in self.ordered_artifacts.all():
        path_list.append(artifact.cache())
    _track_run_input(self, is_run_input)
    return path_list


# docstring handled through attach_func_to_class_method
def load(
    self,
    join: Literal["inner", "outer"] = "outer",
    is_run_input: bool | None = None,
    **kwargs,
) -> Any:
    # cannot call _track_run_input here, see comment further down
    all_artifacts = self.ordered_artifacts.all()
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
def delete(self, permanent: bool | None = None) -> None:
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
def save(self, using: str | None = None) -> Collection:
    if self.meta_artifact is not None:
        self.meta_artifact.save()
    # we don't need to save feature sets again
    save_feature_sets(self)
    super(Collection, self).save()
    # we don't allow updating the collection of artifacts
    # if users want to update the set of artifacts, they
    # have to create a new collection
    if hasattr(self, "_artifacts"):
        links = [
            CollectionArtifact(collection_id=self.id, artifact_id=artifact.id)
            for artifact in self._artifacts
        ]
        # the below seems to preserve the order of the list in the
        # auto-incrementing integer primary
        # merely using .artifacts.set(*...) doesn't achieve this
        # we need ignore_conflicts=True so that this won't error if links already exist
        CollectionArtifact.objects.bulk_create(links, ignore_conflicts=True)
    save_feature_set_links(self)
    if using is not None:
        logger.warning("using argument is ignored")
    return self


# docstring handled through attach_func_to_class_method
def restore(self) -> None:
    self.visibility = VisibilityChoice.default.value
    self.save()


@property  # type: ignore
@doc_args(Collection.ordered_artifacts.__doc__)
def ordered_artifacts(self) -> QuerySet:
    """{}"""  # noqa: D415
    return self.artifacts.order_by("links_collection__id")


@property  # type: ignore
@doc_args(Collection.data_artifact.__doc__)
def data_artifact(self) -> Artifact | None:
    """{}"""  # noqa: D415
    return self.artifacts.first()


METHOD_NAMES = [
    "__init__",
    "mapped",
    "cache",
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

Collection.ordered_artifacts = ordered_artifacts
Collection.data_artifact = data_artifact
