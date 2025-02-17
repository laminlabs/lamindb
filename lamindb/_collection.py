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

from ._parents import view_lineage
from ._record import _get_record_kwargs, init_self_from_db, update_attributes
from ._utils import attach_func_to_class_method
from .core._data import (
    _track_run_input,
    describe,
    get_run,
    save_schema_links,
    save_staged_feature_sets,
)
from .core._mapped_collection import MappedCollection
from .core.storage._pyarrow_dataset import _is_pyarrow_dataset, _open_pyarrow_dataset
from .core.versioning import process_revises
from .errors import FieldValidationError
from .models import (
    Artifact,
    Collection,
    CollectionArtifact,
    Run,
    Schema,
)

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pyarrow.dataset import Dataset as PyArrowDataset

    from ._query_set import QuerySet
    from .core.storage import UPath


class CollectionFeatureManager:
    """Query features of artifact in collection."""

    def __init__(self, collection: Collection):
        self._collection = collection

    def _get_staged_feature_sets_union(self) -> dict[str, Schema]:
        links_schema_artifact = Artifact.feature_sets.through.objects.filter(
            artifact_id__in=self._collection.artifacts.values_list("id", flat=True)
        )
        feature_sets_by_slots = defaultdict(list)
        for link in links_schema_artifact:
            feature_sets_by_slots[link.slot].append(link.schema_id)
        feature_sets_union = {}
        for slot, schema_ids_slot in feature_sets_by_slots.items():
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
            feature_sets_union[slot] = Schema(features, dtype=schema_1.dtype)
        return feature_sets_union


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
    meta_artifact: Artifact | None = kwargs.pop("meta_artifact", None)
    tmp_key: str | None = kwargs.pop("key", None)
    description: str | None = kwargs.pop("description", None)
    reference: str | None = kwargs.pop("reference", None)
    reference_type: str | None = kwargs.pop("reference_type", None)
    run: Run | None = kwargs.pop("run", None)
    revises: Collection | None = kwargs.pop("revises", None)
    version: str | None = kwargs.pop("version", None)
    _branch_code: int | None = kwargs.pop("_branch_code", 1)
    key: str
    if "name" in kwargs:
        key = kwargs.pop("name")
        warnings.warn(
            f"argument `name` will be removed, please pass {key} to `key` instead",
            FutureWarning,
            stacklevel=2,
        )
    else:
        key = tmp_key
    if not len(kwargs) == 0:
        valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Collection)])
        raise FieldValidationError(
            f"Only {valid_keywords} can be passed, you passed: {kwargs}"
        )
    if revises is None:
        revises = (
            Collection.filter(key=key, is_latest=True).order_by("-created_at").first()
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
            f"returning existing collection with same hash: {existing_collection}; if you intended to query to track this collection as an input, use: ln.Collection.get()"
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
        _skip_validation = revises is not None and key == revises.key
        super(Collection, collection).__init__(  # type: ignore
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
            _skip_validation=_skip_validation,
        )
    collection._artifacts = artifacts
    # register provenance
    if revises is not None:
        _track_run_input(revises, run=run)
    _track_run_input(artifacts, run=run)


# docstring handled through attach_func_to_class_method
def append(self, artifact: Artifact, run: Run | None = None) -> Collection:
    return Collection(  # type: ignore
        self.artifacts.all().list() + [artifact],
        # key is automatically taken from revises.key
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
def open(self, is_run_input: bool | None = None) -> PyArrowDataset:
    if self._state.adding:
        artifacts = self._artifacts
        logger.warning("the collection isn't saved, consider calling `.save()`")
    else:
        artifacts = self.ordered_artifacts.all()
    paths = [artifact.path for artifact in artifacts]
    # this checks that the filesystem is the same for all paths
    # this is a requirement of pyarrow.dataset.dataset
    fs = paths[0].fs
    for path in paths[1:]:
        # this assumes that the filesystems are cached by fsspec
        if path.fs is not fs:
            raise ValueError(
                "The collection has artifacts with different filesystems, this is not supported."
            )
    if not _is_pyarrow_dataset(paths):
        suffixes = {path.suffix for path in paths}
        suffixes_str = ", ".join(suffixes)
        err_msg = "This collection is not compatible with pyarrow.dataset.dataset(), "
        err_msg += (
            f"the artifacts have incompatible file types: {suffixes_str}"
            if len(suffixes) > 1
            else f"the file type {suffixes_str} is not supported by pyarrow."
        )
        raise ValueError(err_msg)
    dataset = _open_pyarrow_dataset(paths)
    # track only if successful
    _track_run_input(self, is_run_input)
    return dataset


# docstring handled through attach_func_to_class_method
def mapped(
    self,
    layers_keys: str | list[str] | None = None,
    obs_keys: str | list[str] | None = None,
    obsm_keys: str | list[str] | None = None,
    obs_filter: dict[str, str | list[str]] | None = None,
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
        logger.warning("the collection isn't saved, consider calling `.save()`")
    else:
        artifacts = self.ordered_artifacts.all()
    for artifact in artifacts:
        if artifact.suffix not in {".h5ad", ".zarr"}:
            logger.warning(f"ignoring artifact with suffix {artifact.suffix}")
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
    save_staged_feature_sets(self)
    super(Collection, self).save()
    # we don't allow updating the collection of artifacts
    # if users want to update the set of artifacts, they
    # have to create a new collection
    if hasattr(self, "_artifacts"):
        links = [
            CollectionArtifact(collection_id=self.id, artifact_id=artifact.id)  # type: ignore
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
    "open",
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

# mypy: ignore-errors
Collection.ordered_artifacts = ordered_artifacts
Collection.data_artifact = data_artifact
Collection.describe = describe
Collection.view_lineage = view_lineage
