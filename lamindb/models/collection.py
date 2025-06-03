from __future__ import annotations

import warnings
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    overload,
)

import anndata as ad
import pandas as pd
from django.db import models
from django.db.models import CASCADE, PROTECT
from lamin_utils import logger
from lamindb_setup.core.hashing import HASH_LENGTH, hash_set

from lamindb.base.fields import (
    CharField,
    ForeignKey,
    OneToOneField,
    TextField,
)

from ..base.ids import base62_20
from ..core._mapped_collection import MappedCollection
from ..core.storage._backed_access import _open_dataframe
from ..errors import FieldValidationError
from ..models._is_versioned import process_revises
from ._is_versioned import IsVersioned
from .artifact import (
    Artifact,
    _populate_subsequent_runs_,
    _track_run_input,
    describe_artifact_collection,
    get_run,
    save_schema_links,
)
from .has_parents import view_lineage
from .run import Run, TracksRun, TracksUpdates
from .sqlrecord import (
    BaseSQLRecord,
    IsLink,
    SQLRecord,
    _get_record_kwargs,
    init_self_from_db,
    update_attributes,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from polars import LazyFrame as PolarsLazyFrame
    from pyarrow.dataset import Dataset as PyArrowDataset

    from ..core.storage import UPath
    from .project import Project, Reference
    from .query_set import QuerySet
    from .transform import Transform
    from .ulabel import ULabel


# below is a draft for the future, see also the tests in test_collection.py
#
# class CollectionFeatureManager:
#     """Query features of artifact in collection."""

#     def __init__(self, collection: Collection):
#         self._collection = collection

#     def _get_staged_feature_sets_union(self) -> dict[str, Schema]:
#         links_schema_artifact = Artifact.feature_sets.through.objects.filter(
#             artifact_id__in=self._collection.artifacts.values_list("id", flat=True)
#         )
#         feature_sets_by_slots = defaultdict(list)
#         for link in links_schema_artifact:
#             feature_sets_by_slots[link.slot].append(link.schema_id)
#         feature_sets_union = {}
#         for slot, schema_ids_slot in feature_sets_by_slots.items():
#             schema_1 = Schema.get(id=schema_ids_slot[0])
#             related_name = schema_1._get_related_name()
#             features_registry = getattr(Schema, related_name).field.model
#             # this way of writing the __in statement turned out to be the fastest
#             # evaluated on a link table with 16M entries connecting 500 feature sets with
#             # 60k genes
#             feature_ids = (
#                 features_registry.schemas.through.objects.filter(
#                     schema_id__in=schema_ids_slot
#                 )
#                 .values(f"{features_registry.__name__.lower()}_id")
#                 .distinct()
#             )
#             features = features_registry.filter(id__in=feature_ids)
#             feature_sets_union[slot] = Schema(features, dtype=schema_1.dtype)
#         return feature_sets_union


def _load_concat_artifacts(
    artifacts: list[Artifact], join: Literal["inner", "outer"] = "outer", **kwargs
) -> pd.DataFrame | ad.AnnData:
    suffixes = {artifact.suffix for artifact in artifacts}
    # Why is that? - Sergei
    if len(suffixes) != 1:
        raise ValueError(
            "Can only load collections where all artifacts have the same suffix"
        )

    # because we're tracking data flow on the collection-level, here, we don't
    # want to track it on the artifact-level
    first_object = artifacts[0].load(is_run_input=False)
    is_dataframe = isinstance(first_object, pd.DataFrame)
    is_anndata = isinstance(first_object, ad.AnnData)
    if not is_dataframe and not is_anndata:
        raise ValueError(f"Unable to concatenate {suffixes.pop()} objects.")

    objects = [first_object]
    artifact_uids = [artifacts[0].uid]
    for artifact in artifacts[1:]:
        objects.append(artifact.load(is_run_input=False))
        artifact_uids.append(artifact.uid)

    if is_dataframe:
        concat_object = pd.concat(objects, join=join, **kwargs)
    elif is_anndata:
        label = kwargs.pop("label", "artifact_uid")
        keys = kwargs.pop("keys", artifact_uids)
        concat_object = ad.concat(objects, join=join, label=label, keys=keys, **kwargs)
    return concat_object


class Collection(SQLRecord, IsVersioned, TracksRun, TracksUpdates):
    """Versioned collections of artifacts.

    Args:
        artifacts: `list[Artifact]` A list of artifacts.
        key: `str` A file-path like key, analogous to the `key` parameter of `Artifact` and `Transform`.
        description: `str | None = None` A description.
        revises: `Collection | None = None` An old version of the collection.
        run: `Run | None = None` The run that creates the collection.
        meta: `Artifact | None = None` An artifact that defines metadata for the collection.
        reference: `str | None = None` A simple reference, e.g. an external ID or a URL.
        reference_type: `str | None = None` A way to indicate to indicate the type of the simple reference `"url"`.

    See Also:
        :class:`~lamindb.Artifact`

    Examples:

        Create a collection from a list of :class:`~lamindb.Artifact` objects:

        >>> collection = ln.Collection([artifact1, artifact2], key="my_project/my_collection")

        Create a collection that groups a data & a metadata artifact (e.g., here :doc:`docs:rxrx`):

        >>> collection = ln.Collection(data_artifact, key="my_project/my_collection", meta=metadata_artifact)

    """

    class Meta(SQLRecord.Meta, IsVersioned.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _len_full_uid: int = 20
    _len_stem_uid: int = 16
    _name_field: str = "key"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False,
        unique=True,
        db_index=True,
        max_length=_len_full_uid,
        default=base62_20,
    )
    """Universal id, valid across DB instances."""
    key: str = CharField(db_index=True)
    """Name or path-like key."""
    # below is the only case in which we use a TextField
    # for description; we do so because users had descriptions exceeding 255 chars
    # in their instances
    description: str | None = TextField(null=True, db_index=True)
    """A description or title."""
    hash: str | None = CharField(
        max_length=HASH_LENGTH, db_index=True, null=True, unique=True
    )
    """Hash of collection content."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A reference like URL or external ID."""
    # also for reference_type here, we allow an extra long max_length
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of reference, e.g., cellxgene Census collection_id."""
    ulabels: ULabel = models.ManyToManyField(
        "ULabel", through="CollectionULabel", related_name="collections"
    )
    """ULabels sampled in the collection (see :class:`~lamindb.Feature`)."""
    run: Run | None = ForeignKey(
        Run, PROTECT, related_name="output_collections", null=True, default=None
    )
    """:class:`~lamindb.Run` that created the `collection`."""
    input_of_runs: Run = models.ManyToManyField(Run, related_name="input_collections")
    """Runs that use this collection as an input."""
    _subsequent_runs: Run = models.ManyToManyField(
        "Run",
        related_name="_recreated_collections",
        db_table="lamindb_collection__previous_runs",  # legacy name, change in lamindb v2
    )
    """Runs that re-created the record after initial creation."""
    artifacts: Artifact = models.ManyToManyField(
        "Artifact", related_name="collections", through="CollectionArtifact"
    )
    """Artifacts in collection."""
    meta_artifact: Artifact | None = OneToOneField(
        "Artifact",
        PROTECT,
        null=True,
        unique=True,
        related_name="_meta_of_collection",
    )
    """An artifact that stores metadata that indexes a collection.

    It has a 1:1 correspondence with an artifact. If needed, you can access the
    collection from the artifact via a private field:
    `artifact._meta_of_collection`.
    """
    _actions: Artifact = models.ManyToManyField(Artifact, related_name="+")
    """Actions to attach for the UI."""
    projects: Project
    """Linked projects."""
    references: Reference
    """Linked references."""

    @overload
    def __init__(
        self,
        artifacts: list[Artifact],
        key: str,
        description: str | None = None,
        meta: Any | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
        run: Run | None = None,
        revises: Collection | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
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
        branch_id: int | None = kwargs.pop("branch_id", 1)
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
            valid_keywords = ", ".join(
                [val[0] for val in _get_record_kwargs(Collection)]
            )
            raise FieldValidationError(
                f"Only {valid_keywords} can be passed, you passed: {kwargs}"
            )
        if revises is None:
            revises = (
                Collection.filter(key=key, is_latest=True)
                .order_by("-created_at")
                .first()
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
            if run is not None:
                existing_collection._populate_subsequent_runs(run)
            init_self_from_db(self, existing_collection)
            update_attributes(self, {"description": description, "key": key})
        else:
            _skip_validation = revises is not None and key == revises.key
            super().__init__(  # type: ignore
                uid=provisional_uid,
                key=key,
                description=description,
                reference=reference,
                reference_type=reference_type,
                meta_artifact=meta_artifact,
                hash=hash,
                run=run,
                version=version,
                branch_id=branch_id,
                revises=revises,
                _skip_validation=_skip_validation,
            )
        self._artifacts = artifacts
        if revises is not None:
            _track_run_input(revises, run=run)
        _track_run_input(artifacts, run=run)

    @classmethod
    def get(
        cls,
        idlike: int | str | None = None,
        *,
        is_run_input: bool | Run = False,
        **expressions,
    ) -> Artifact:
        """Get a single collection.

        Args:
            idlike: Either a uid stub, uid or an integer id.
            is_run_input: Whether to track this collection as run input.
            expressions: Fields and values passed as Django query expressions.

        Raises:
            :exc:`docs:lamindb.errors.DoesNotExist`: In case no matching record is found.

        See Also:
            - Method in `SQLRecord` base class: :meth:`~lamindb.models.SQLRecord.get`

        Examples:

            ::

                collection = ln.Collection.get("okxPW6GIKBfRBE3B0000")
                collection = ln.Collection.get(key="scrna/collection1")
        """
        from .query_set import QuerySet

        return QuerySet(model=cls).get(idlike, is_run_input=is_run_input, **expressions)

    def append(self, artifact: Artifact, run: Run | None = None) -> Collection:
        """Append an artifact to the collection.

        This does not modify the original collection in-place, but returns a new version
        of the original collection with the appended artifact.

        Args:
            artifact: An artifact to add to the collection.
            run: The run that creates the new version of the collection.

        Examples:

            ::

                collection_v1 = ln.Collection(artifact, key="My collection").save()
                collection_v2 = collection.append(another_artifact)  # returns a new version of the collection
                collection_v2.save()  # save the new version

        """
        return Collection(  # type: ignore
            self.artifacts.all().list() + [artifact],
            # key is automatically derived from revises.key
            description=self.description,
            revises=self,
            run=run,
        )

    def open(
        self,
        engine: Literal["pyarrow", "polars"] = "pyarrow",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> PyArrowDataset | Iterator[PolarsLazyFrame]:
        """Open a dataset for streaming.

        Works for `pyarrow` and `polars` compatible formats
        (`.parquet`, `.csv`, `.ipc` etc. files or directories with such files).

        Args:
            engine: Which module to use for lazy loading of a dataframe
                from `pyarrow` or `polars` compatible formats.
            is_run_input: Whether to track this artifact as run input.
            **kwargs: Keyword arguments for `pyarrow.dataset.dataset` or `polars.scan_*` functions.

        Notes:
            For more info, see guide: :doc:`/arrays`.
        """
        if self._state.adding:
            artifacts = self._artifacts
            logger.warning("the collection isn't saved, consider calling `.save()`")
        else:
            artifacts = self.ordered_artifacts.all()
        paths = [artifact.path for artifact in artifacts]

        dataframe = _open_dataframe(paths, engine=engine, **kwargs)
        # track only if successful
        # is it really needed if tracking is done in self.ordered_artifacts.all()? - Sergei
        _track_run_input(self, is_run_input)
        return dataframe

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
        """Return a map-style dataset.

        Returns a `pytorch map-style dataset
        <https://pytorch.org/docs/stable/data.html#map-style-datasets>`__ by
        virtually concatenating `AnnData` arrays.

        By default (`stream=False`) `AnnData` arrays are moved into a local
        cache first.

        `__getitem__` of the `MappedCollection` object takes a single integer index
        and returns a dictionary with the observation data sample for this index from
        the `AnnData` objects in the collection. The dictionary has keys for `layers_keys`
        (`.X` is in `"X"`), `obs_keys`, `obsm_keys` (under `f"obsm_{key}"`) and also `"_store_idx"`
        for the index of the `AnnData` object containing this observation sample.

        .. note::

            For a guide, see :doc:`docs:scrna-mappedcollection`.

            This method currently only works for collections or query sets of `AnnData` artifacts.

        Args:
            layers_keys: Keys from the ``.layers`` slot. ``layers_keys=None`` or ``"X"`` in the list
                retrieves ``.X``.
            obs_keys: Keys from the ``.obs`` slots.
            obsm_keys: Keys from the ``.obsm`` slots.
            obs_filter: Select only observations with these values for the given obs columns.
                Should be a dictionary with obs column names as keys
                and filtering values (a string or a list of strings) as values.
            join: `"inner"` or `"outer"` virtual joins. If ``None`` is passed,
                does not join.
            encode_labels: Encode labels into integers.
                Can be a list with elements from ``obs_keys``.
            unknown_label: Encode this label to -1.
                Can be a dictionary with keys from ``obs_keys`` if ``encode_labels=True``
                or from ``encode_labels`` if it is a list.
            cache_categories: Enable caching categories of ``obs_keys`` for faster access.
            parallel: Enable sampling with multiple processes.
            dtype: Convert numpy arrays from ``.X``, ``.layers`` and ``.obsm``
            stream: Whether to stream data from the array backend.
            is_run_input: Whether to track this collection as run input.

        Examples:
            >>> import lamindb as ln
            >>> from torch.utils.data import DataLoader
            >>> ds = ln.Collection.get(description="my collection")
            >>> mapped = collection.mapped(obs_keys=["cell_type", "batch"])
            >>> dl = DataLoader(mapped, batch_size=128, shuffle=True)
            >>> # also works for query sets of artifacts, '...' represents some filtering condition
            >>> # additional filtering on artifacts of the collection
            >>> mapped = collection.artifacts.all().filter(...).order_by("-created_at").mapped()
            >>> # or directly from a query set of artifacts
            >>> mapped = ln.Artifact.filter(..., otype="AnnData").order_by("-created_at").mapped()
        """
        path_list = []
        if self._state.adding:
            artifacts = self._artifacts
            logger.warning("the collection isn't saved, consider calling `.save()`")
        else:
            artifacts = self.ordered_artifacts.all()
        for artifact in artifacts:
            if ".h5ad" not in artifact.suffix and ".zarr" not in artifact.suffix:
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
        # is it really needed if tracking is done in self.ordered_artifacts.all()? - Sergei
        _track_run_input(self, is_run_input)
        return ds

    def cache(self, is_run_input: bool | None = None) -> list[UPath]:
        """Download cloud artifacts in collection to local cache.

        Follows synching logic: only caches outdated artifacts.

        Returns paths to locally cached on-disk artifacts.

        Args:
            is_run_input: Whether to track this collection as run input.
        """
        path_list = []
        for artifact in self.ordered_artifacts.all():
            path_list.append(artifact.cache())
        # is it really needed if tracking is done in self.ordered_artifacts.all()? - Sergei
        _track_run_input(self, is_run_input)
        return path_list

    def load(
        self,
        join: Literal["inner", "outer"] = "outer",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> pd.DataFrame | ad.AnnData:
        """Cache and load to memory.

        Returns an in-memory concatenated `DataFrame` or `AnnData` object.
        """
        # cannot call _track_run_input here, see comment further down
        artifacts = self.ordered_artifacts.all()
        concat_object = _load_concat_artifacts(artifacts, join, **kwargs)
        # only call it here because there might be errors during load or concat
        # is it really needed if tracking is done in self.ordered_artifacts.all()? - Sergei
        _track_run_input(self, is_run_input)
        return concat_object

    def delete(self, permanent: bool | None = None) -> None:
        """Delete collection.

        Args:
            permanent: Whether to permanently delete the collection record (skips trash).

        Examples:

            For any `Collection` object `collection`, call:

            >>> collection.delete()
        """
        # change branch_id to trash
        trash_branch_id = -1
        if self.branch_id > trash_branch_id and permanent is not True:
            self.branch_id = trash_branch_id
            self.save()
            logger.warning(f"moved collection to trash (branch_id = {trash_branch_id})")
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
            super().delete()

    def save(self, using: str | None = None) -> Collection:
        """Save the collection and underlying artifacts to database & storage.

        Args:
            using: The database to which you want to save.

        Examples:
            >>> collection = ln.Collection("./myfile.csv", name="myfile")
        """
        if self.meta_artifact is not None:
            self.meta_artifact.save()
        super().save()
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

    def restore(self) -> None:
        """Restore collection record from trash.

        Examples:

            For any `Collection` object `collection`, call:

            >>> collection.restore()
        """
        self.branch_id = 1
        self.save()

    @property
    def transform(self) -> Transform | None:
        """Transform whose run created the collection."""
        return self.run.transform if self.run is not None else None

    @property
    def name(self) -> str:
        """Name of the collection.

        Splits `key` on `/` and returns the last element.
        """
        return self.key.split("/")[-1]

    @property
    def ordered_artifacts(self) -> QuerySet:
        """Ordered `QuerySet` of `.artifacts`.

        Accessing the many-to-many field `collection.artifacts` directly gives
        you non-deterministic order.

        Using the property `.ordered_artifacts` allows to iterate through a set
        that's ordered in the order of creation.
        """
        return self.artifacts.order_by("links_collection__id")

    @property
    def data_artifact(self) -> Artifact | None:
        """Access to a single data artifact.

        If the collection has a single data & metadata artifact, this allows access via::

           collection.data_artifact  # first & only element of collection.artifacts
           collection.meta_artifact  # metadata

        """
        return self.artifacts.first()

    def describe(self) -> None:
        """Describe relations of record.

        Examples:
            >>> artifact.describe()
        """
        return describe_artifact_collection(self)

    def _populate_subsequent_runs(self, run: Run) -> None:
        _populate_subsequent_runs_(self, run)


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


class CollectionArtifact(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_artifact"
    )
    artifact: Artifact = ForeignKey(Artifact, PROTECT, related_name="links_collection")

    class Meta:
        unique_together = ("collection", "artifact")


# mypy: ignore-errors
Collection.view_lineage = view_lineage
