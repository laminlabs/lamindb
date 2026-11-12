from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal

from django.db import transaction
from lamindb_setup import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core.upath import LocalPathClasses

from lamindb.errors import IntegrityError

from ._is_versioned import create_uid

if TYPE_CHECKING:
    from polars import LazyFrame
    from upath import UPath

    from .artifact import Artifact
    from .query_set import QuerySet
    from .schema import Schema


class StreamedLazyArtifact:
    """Managed lazy writes with in-memory UID allocation, not reservations."""

    def __init__(
        self,
        *,
        key: str | None,
        frame: LazyFrame | None = None,
        schema: Schema | None = None,
        parquet_kwargs: dict[str, Any] | None = None,
        on_validation_error: Literal["raise", "keep"] = "raise",
        **kwargs,
    ):
        from .artifact import Artifact
        from .sqlrecord import check_key, get_current_branch

        if key is None and frame is None:
            raise ValueError(
                "Caller-managed streaming requires an explicit virtual key."
            )
        if key is not None:
            check_key(key)
        if kwargs.pop("key_is_virtual", True) is not True:
            raise ValueError("Streamed lazy artifacts require virtual keys.")
        allowed = {
            "description",
            "run",
            "revises",
            "storage",
            "branch",
            "version_tag",
            "features",
            "otype",
            "kind",
        }
        unsupported = set(kwargs) - allowed
        if unsupported:
            raise ValueError(f"Unsupported streaming arguments: {sorted(unsupported)}")
        storage = kwargs.pop("storage", None)
        if storage is None:
            storage = setup_settings.storage.record
            from lamindb import context

            space = context.space or setup_settings.space
            if space is not None and space.id != storage.space_id:
                raise ValueError(
                    "Pass storage in the active space explicitly for streaming."
                )
        if storage._state.db not in {None, "default"}:
            raise ValueError("Streaming supports the current database only.")
        if storage.instance_uid != setup_settings.instance.uid:
            raise ValueError("Streaming requires storage managed by this instance.")
        self._root = storage.path
        if not isinstance(self._root, LocalPathClasses) and self._root.protocol != "s3":
            raise ValueError("Streaming currently supports S3 and local storage only.")
        self._storage = storage
        self._branch = kwargs.pop("branch", None) or get_current_branch()
        self._key = key
        self._kwargs = kwargs
        self._result: Artifact | None = None
        self._metadata_ready = False
        self._cleanup_pending = False
        self._discarded = False
        self._frame = frame
        self._validated_frame = None
        self._on_validation_error = on_validation_error
        self._sink_options = dict(parquet_kwargs or {})
        self._sink_complete = False
        self._sink_failed = False
        if self._sink_options.get("lazy") or "path" in self._sink_options:
            raise ValueError("The sink must execute at the managed path.")
        if frame is not None and schema is not None:
            from lamindb.curators._polars_streaming import ValidatedLazyFrame

            if schema.suffix not in {None, ".parquet"}:
                raise ValueError("The schema suffix must match the Parquet sink.")
            self._validated_frame = ValidatedLazyFrame(
                frame,
                schema,
                features=kwargs.get("features"),
                on_validation_error=on_validation_error,
            )
            self._frame = self._validated_frame.frame
        revises = kwargs.get("revises")
        if revises is not None and (
            not isinstance(revises, Artifact)
            or revises.pk is None
            or revises.storage_id != storage.id
            or revises.key != key
            or revises._state.db not in {None, "default"}
        ):
            raise ValueError(
                "revises must be a saved artifact in this database "
                "with the same storage and key."
            )
        if transaction.get_connection().in_atomic_block:
            raise ValueError(
                "Create streamed lazy artifacts outside an existing database transaction."
            )
        self._allocate_path()

    @property
    def validation_error(self) -> str | None:
        if self._validated_frame is None:
            return None
        return self._validated_frame.error

    def _head(self) -> Artifact | None:
        records = self._artifact_records()
        if self._key is None:
            revises = self._kwargs.get("revises")
            if revises is None:
                return None
            records = records.filter(uid__startswith=revises.stem_uid)
        heads = records.filter(is_latest=True)
        head = (
            heads.filter(branch=self._branch).order_by("-created_at", "-id").first()
            or heads.order_by("-created_at", "-id").first()
        )
        if head is None and records.exists():
            raise IntegrityError(
                f"Artifacts exist for {self._key!r}, but none is a latest version."
            )
        return head

    def _artifact_records(self) -> QuerySet:
        from .artifact import Artifact

        return Artifact.objects.filter(storage=self._storage, key=self._key).exclude(
            branch_id=-1
        )

    def _make_artifact(
        self,
        uid: str,
        path: UPath,
        revises: Artifact | None,
        stats: dict[str, Any],
    ) -> Artifact:
        from .artifact import Artifact

        options = {k: v for k, v in self._kwargs.items() if k != "revises"}
        options.setdefault("otype", "DataFrame")
        options.setdefault("kind", "dataset")
        artifact = Artifact(
            path,
            key=self._key,
            storage=self._storage,
            space=self._storage.space,
            branch=self._branch,
            revises=revises,
            overwrite_versions=False,
            _is_internal_call=True,
            skip_check_exists=True,
            skip_hash_lookup=True,
            _stats=stats,
            **options,
        )
        if artifact.uid != uid:
            raise IntegrityError("Managed destination does not match artifact UID.")
        artifact._refresh_revises_if_stale = False
        return artifact

    def _allocate_path(self) -> None:
        from lamindb.core.storage.paths import auto_storage_key_from_artifact_uid

        head = self._head()
        explicit_revises = self._kwargs.get("revises")
        if (
            explicit_revises is not None
            and head is not None
            and explicit_revises.stem_uid != head.stem_uid
        ):
            raise IntegrityError(
                "revises does not belong to the current key's version family."
            )
        revises = head or explicit_revises
        uid = create_uid(revises=revises, version_tag=self._kwargs.get("version_tag"))
        path = self._root / auto_storage_key_from_artifact_uid(
            uid, ".parquet", overwrite_versions=False
        )
        if path.exists():
            raise IntegrityError(
                f"Streaming destination already exists at {path}. "
                "Clean up its unfinished writer before starting another write."
            )
        self._baseline_uid = revises.uid if revises is not None else None
        self._artifact = self._make_artifact(uid, path, revises, {})
        self._uid = uid
        self.path: UPath = path
        if isinstance(path, LocalPathClasses):
            path.parent.mkdir(parents=True, exist_ok=True)

    def _read_metadata(self) -> None:
        import pyarrow.parquet as pq

        from .artifact import get_stat_file_cloud, hash_file

        path = self.path
        if path.is_dir():
            raise ValueError("Streaming currently supports single-file Parquet only.")
        if isinstance(path, LocalPathClasses):
            size, hash_value, hash_type = hash_file(path)
        else:
            path.fs.invalidate_cache(str(path))
            size, hash_value, hash_type = get_stat_file_cloud(
                path.stat().as_info(), protocol=path.protocol
            )
        if hash_value is None:
            raise ValueError(
                f"Storage did not provide a hash for {path}; upload retained."
            )
        with path.open("rb") as file:
            n_observations = pq.read_metadata(file).num_rows
        self._artifact.size = size
        self._artifact.hash = hash_value
        self._artifact._hash_type = hash_type
        self._artifact.n_observations = n_observations
        self._metadata_ready = True

    def _sink(self) -> None:
        if self.path.exists():
            raise IntegrityError(
                f"Streaming destination already exists at {self.path}; "
                "refusing to overwrite it. Writers must be serialized."
            )
        options = dict(self._sink_options)
        if (
            self.path.protocol == "s3"
            and "storage_options" not in options
            and "credential_provider" not in options
        ):
            from lamindb.core.storage._polars_lazy_df import _polars_options

            options.update(_polars_options(self.path))
        try:
            self._frame.sink_parquet(self.path.as_posix(), **options)
        except Exception as sink_error:
            self._sink_failed = True
            if (
                self._on_validation_error == "raise"
                and self.validation_error is not None
            ):
                logger.error(
                    f"Validation failed at {self.path}; no artifact was published. "
                    "Removing partial output."
                )
                try:
                    self.discard()
                except Exception as cleanup_error:
                    logger.error(
                        f"Validation cleanup failed at {self.path}: {cleanup_error}. "
                        "Use discard() to retry cleanup."
                    )
                    raise sink_error from cleanup_error
            else:
                logger.error(
                    f"Sinking failed at {self.path}; no artifact was published. "
                    "Use discard() to remove any partial object."
                )
            raise
        self._sink_complete = True

    def save(self) -> Artifact:
        """Deduplicate and publish an already-written file without renaming it."""
        if self._discarded:
            raise ValueError("This streamed write has been discarded.")
        if self._result is not None:
            if self._cleanup_pending:
                self._cleanup()
                self._cleanup_pending = False
            return self._result
        if transaction.get_connection().in_atomic_block:
            raise ValueError(
                "Streaming save must run outside an existing database transaction."
            )
        if self._sink_failed:
            raise ValueError("The sink failed; discard this write and start a new one.")
        if self._artifact.__class__.objects.filter(uid=self._uid).exists():
            raise IntegrityError(
                f"The streaming UID {self._uid} was published by another writer; "
                "refusing to finalize or remove its object. Writers must be serialized."
            )
        if self._frame is not None and not self._sink_complete:
            self._sink()
        if not self._metadata_ready:
            self._read_metadata()
        with transaction.atomic():
            match = (
                self._artifact_records()
                .filter(hash=self._artifact.hash)
                .order_by("-created_at", "-id")
                .first()
            )
            if match is not None:
                result = match
                duplicate = True
            else:
                head = self._head()
                if (head.uid if head is not None else None) != self._baseline_uid:
                    raise IntegrityError(
                        f"A different version of {self._key!r} was published during "
                        f"streaming. Upload retained at {self.path}; "
                        "discard it and start a new write to revise the current head."
                    )
                artifact = self._make_artifact(
                    self._uid,
                    self.path,
                    head,
                    {
                        "size": self._artifact.size,
                        "hash": self._artifact.hash,
                        "hash_type": self._artifact._hash_type,
                    },
                )
                artifact.n_observations = self._artifact.n_observations
                if self._validated_frame is not None:
                    curator = self._validated_frame.curator
                    if self.validation_error is None:
                        curator._artifact = artifact
                        artifact.schema = curator._schema
                        artifact._curator = curator
                    else:
                        artifact.extra_data = {
                            **(artifact.extra_data or {}),
                            "validation": {
                                "status": "failed",
                                "schema_uid": curator._schema.uid,
                                "error": self.validation_error,
                            },
                        }
                result = artifact.save()
                if result.uid != self._uid:
                    raise IntegrityError(
                        "Finalization changed the streaming artifact UID."
                    )
                duplicate = False
        self._result = result
        if duplicate:
            self._cleanup_pending = True
            self._cleanup()
            self._cleanup_pending = False
        return result

    def _cleanup(self) -> None:
        if self._artifact.__class__.objects.filter(uid=self._uid).exists():
            raise ValueError(
                "Cannot discard a published artifact; use Artifact.delete()."
            )
        try:
            self.path.unlink(missing_ok=True)
        except Exception:
            logger.error(
                f"Could not remove redundant or abandoned streaming object at "
                f"{self.path}; use discard() to retry cleanup."
            )
            raise

    def discard(self) -> None:
        """Remove an unfinished write, never an object with a published UID."""
        if self.path.is_dir():
            raise ValueError("Folder cleanup is not supported by this prototype.")
        self._cleanup()
        self._discarded = True
