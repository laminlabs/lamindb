from __future__ import annotations

import warnings
from collections import defaultdict
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
)

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.hashing import hash_set

from lamindb.models import (
    Collection,
    CollectionArtifact,
    Schema,
)

from ._parents import view_lineage
from ._record import init_self_from_db, update_attributes
from ._utils import attach_func_to_class_method
from .core._data import (
    _track_run_input,
    describe,
    get_run,
    save__schemas_m2m,
    save_schema_links,
)
from .core._mapped_collection import MappedCollection
from .core._settings import settings
from .core.versioning import process_revises
from .models import Artifact, Run

if TYPE_CHECKING:
    from collections.abc import Iterable

    from lamindb.core.storage import UPath

    from ._query_set import QuerySet


class CollectionFeatureManager:
    """Query features of artifact in collection."""

    def __init__(self, collection: Collection):
        self._collection = collection

    def _get__schemas_m2m_union(self) -> dict[str, Schema]:
        links_schema_artifact = Artifact._schemas_m2m.through.objects.filter(
            artifact_id__in=self._collection.artifacts.values_list("id", flat=True)
        )
        _schemas_m2m_by_slots = defaultdict(list)
        for link in links_schema_artifact:
            _schemas_m2m_by_slots[link.slot].append(link.schema_id)
        _schemas_m2m_union = {}
        for slot, schema_ids_slot in _schemas_m2m_by_slots.items():
            schema_1 = Schema.get(id=schema_ids_slot[0])
            related_name = schema_1._get_related_name()
            features_registry = getattr(Schema, related_name).field.model
            # this way of writing the __in statement turned out to be the fastest
            # evaluated on a link table with 16M entries connecting 500 feature sets with
            # 60k genes
            feature_ids = (
                features_registry.schemas.through.objects.filter(
                    schema_id__in=schema_ids_slot
                )
                .values(f"{features_registry.__name__.lower()}_id")
                .distinct()
            )
            features = features_registry.filter(id__in=feature_ids)
            _schemas_m2m_union[slot] = Schema(features, dtype=schema_1.dtype)
        return _schemas_m2m_union


def __init__(
    collection: Collection,
    *args,
    **kwargs,
):
    collection.features = CollectionFeatureManager(collection)
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
    key: str | None = kwargs.pop("key") if "key" in kwargs else None
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
    _branch_code: int | None = (
        kwargs.pop("_branch_code") if "_branch_code" in kwargs else 1
    )
    if "name" in kwargs:
        key = kwargs.pop("name")
        warnings.warn(
            f"argument `name` will be removed, please pass {key} to `key` instead",
            FutureWarning,
            stacklevel=2,
        )
    if not len(kwargs) == 0:
        raise ValueError(
            f"Only artifacts, key, run, description, reference, reference_type can be passed, you passed: {kwargs}"
        )
    provisional_uid, version, key, description, revises = process_revises(
        revises, version, key, description, Collection
    )
    run = get_run(run)
    if isinstance(artifacts, Artifact):
        artifacts = [artifacts]
    else:
        if not hasattr(artifacts, "__getitem__"):
            raise ValueError("Artifact or list[Artifact] is allowed.")
        assert isinstance(artifacts[0], Artifact)  # type: ignore  # noqa: S101
    hash = from_artifacts(artifacts)  # type: ignore
    if meta_artifact is not None:
        if not isinstance(meta_artifact, Artifact):
            raise ValueError("meta_artifact has to be an Artifact")
        if isinstance(meta_artifact, Artifact):
            if meta_artifact._state.adding:
                raise ValueError(
                    "Save meta_artifact artifact before creating collection!"
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
        # update the run of the existing collection
        if run is not None:
            # save the information that this collection was previously produced
            # by another run
            # note: same logic exists for _output_artifacts_with_later_updates
            if existing_collection.run is not None and existing_collection.run != run:
                existing_collection.run._output_collections_with_later_updates.add(
                    existing_collection
                )
            # update the run of the collection with the latest run
            existing_collection.run = run
        init_self_from_db(collection, existing_collection)
        update_attributes(collection, {"description": description, "key": key})
    else:
        kwargs = {}
        search_names_setting = settings.creation.search_names
        if revises is not None and key == revises.key:
            settings.creation.search_names = False
        super(Collection, collection).__init__(
            uid=provisional_uid,
            key=key,
            description=description,
            reference=reference,
            reference_type=reference_type,
            meta_artifact=meta_artifact,
            hash=hash,
            run=run,
            version=version,
            _branch_code=_branch_code,
            revises=revises,
            **kwargs,
        )
        settings.creation.search_names = search_names_setting
    collection._artifacts = artifacts
    # register provenance
    if revises is not None:
        _track_run_input(revises, run=run)
    _track_run_input(artifacts, run=run)


# docstring handled through attach_func_to_class_method
def append(self, artifact: Artifact, run: Run | None = None) -> Collection:
    return Collection(
        self.artifacts.all().list() + [artifact],
        description=self.description,
        revises=self,
        run=run,
    )


# internal function, not exposed to user
def from_artifacts(artifacts: Iterable[Artifact]) -> tuple[str, dict[str, str]]:
    # assert all artifacts are already saved
    saved = not any(artifact._state.adding for artifact in artifacts)
    if not saved:
        raise ValueError("Not all artifacts are yet saved, please save them")
    # validate consistency of hashes - we do not allow duplicate hashes
    hashes = [artifact.hash for artifact in artifacts if artifact.hash is not None]
    hashes_set = set(hashes)
    if len(hashes) != len(hashes_set):
        seen = set()
        non_unique = [x for x in hashes if x in seen or seen.add(x)]  # type: ignore
        raise ValueError(
            "Please pass artifacts with distinct hashes: these ones are non-unique"
            f" {non_unique}"
        )
    hash = hash_set(hashes_set)
    return hash


# docstring handled through attach_func_to_class_method
def mapped(
    self,
    layers_keys: str | list[str] | None = None,
    obs_keys: str | list[str] | None = None,
    obsm_keys: str | list[str] | None = None,
    obs_filter: tuple[str, str | tuple[str, ...]] | None = None,
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
    if self._state.adding:
        artifacts = self._artifacts
        logger.warning("The collection isn't saved, consider calling `.save()`")
    else:
        artifacts = self.ordered_artifacts.all()
    for artifact in artifacts:
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
        obs_filter,
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
    # change _branch_code to trash
    trash__branch_code = -1
    if self._branch_code > trash__branch_code and permanent is not True:
        self._branch_code = trash__branch_code
        self.save()
        logger.warning(
            f"moved collection to trash (_branch_code = {trash__branch_code})"
        )
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
    save__schemas_m2m(self)
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
    save_schema_links(self)
    if using is not None:
        logger.warning("using argument is ignored")
    return self


# docstring handled through attach_func_to_class_method
def restore(self) -> None:
    self._branch_code = 1
    self.save()


@property  # type: ignore
@doc_args(Collection.ordered_artifacts.__doc__)
def ordered_artifacts(self) -> QuerySet:
    """{}"""  # noqa: D415
    # tracking is done via QueryManager (_query_manager.py)
    return self.artifacts.order_by("links_collection__id")


@property  # type: ignore
@doc_args(Collection.data_artifact.__doc__)
def data_artifact(self) -> Artifact | None:
    """{}"""  # noqa: D415
    return self.artifacts.first()


METHOD_NAMES = [
    "__init__",
    "append",
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
Collection.describe = describe
Collection.view_lineage = view_lineage
