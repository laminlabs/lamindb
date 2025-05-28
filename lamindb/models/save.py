# ruff: noqa: TC004
from __future__ import annotations

import os
import shutil
import traceback
from collections import defaultdict
from datetime import datetime
from typing import TYPE_CHECKING

from django.db import transaction
from django.utils.functional import partition
from lamin_utils import logger
from lamindb_setup.core.upath import LocalPathClasses, UPath

from ..core._settings import settings
from ..core.storage.paths import (
    _cache_key_from_artifact_storage,
    attempt_accessing_path,
    auto_storage_key_from_artifact,
    delete_storage_using_key,
    store_file_or_folder,
)
from .sqlrecord import SQLRecord

if TYPE_CHECKING:
    from collections.abc import Iterable

    from .artifact import Artifact


def save(records: Iterable[SQLRecord], ignore_conflicts: bool | None = False) -> None:
    """Bulk save records.

    Note:

        This is a much faster than saving records using ``record.save()``.

    Warning:

        Bulk saving neither automatically creates related records nor updates
        existing records! Use ``record.save()`` for these use cases.

    Args:
        records: Multiple :class:`~lamindb.models.SQLRecord` objects.
        ignore_conflicts: If ``True``, do not error if some records violate a
           unique or another constraint. However, it won't inplace update the id
           fields of records. If you need records with ids, you need to query
           them from the database.

    Examples:

        Save a list of records:

        >>> labels = [ln.ULabel(f"Label {i}") for i in range(10)]
        >>> ln.save(projects)

        For a single record, use ``record.save()``:

        >>> transform = ln.Transform(key="My pipeline")
        >>> transform.save()

        Update a single existing record:

        >>> transform = ln.Transform.get("0Cb86EZj")
        >>> transform.description = "New description"
        >>> transform.save()

    """
    from .artifact import Artifact

    if isinstance(records, SQLRecord):
        raise ValueError("Please use record.save() if saving a single record.")

    # previously, this was all set based,
    # but models without primary keys aren't hashable
    # we distinguish between artifacts and non-artifacts
    # for artifacts, we want to bulk-upload rather than upload one-by-one
    non_artifacts, artifacts = partition(lambda r: isinstance(r, Artifact), records)
    if non_artifacts:
        non_artifacts_old, non_artifacts_new = partition(
            lambda r: r._state.adding or r.pk is None, non_artifacts
        )
        bulk_create(non_artifacts_new, ignore_conflicts=ignore_conflicts)
        if non_artifacts_old:
            bulk_update(non_artifacts_old)
        non_artifacts_with_parents = [
            r for r in non_artifacts_new if hasattr(r, "_parents")
        ]
        if len(non_artifacts_with_parents) > 0:
            # this can only happen within bionty right now!!
            # we might extend to core lamindb later
            from bionty.core import add_ontology

            add_ontology(non_artifacts_with_parents)

    if artifacts:
        with transaction.atomic():
            for record in artifacts:
                record._save_skip_storage()
        using_key = settings._using_key
        store_artifacts(artifacts, using_key=using_key)

    # this function returns None as potentially 10k records might be saved
    # refreshing all of them from the DB would mean a severe performance penalty
    # 2nd reason: consistency with Django Model.save(), which also returns None
    return None


def bulk_create(records: Iterable[SQLRecord], ignore_conflicts: bool | None = False):
    records_by_orm = defaultdict(list)
    for record in records:
        records_by_orm[record.__class__].append(record)
    for registry, records in records_by_orm.items():
        registry.objects.bulk_create(records, ignore_conflicts=ignore_conflicts)
        # records[:] = created  # In-place list update; does not seem to be necessary


def bulk_update(records: Iterable[SQLRecord], ignore_conflicts: bool | None = False):
    records_by_orm = defaultdict(list)
    for record in records:
        records_by_orm[record.__class__].append(record)
    for registry, records in records_by_orm.items():
        field_names = [
            field.name
            for field in registry._meta.fields
            if (field.name != "created_at" and field.name != "id")
        ]
        registry.objects.bulk_update(records, field_names)


# This is also used within Artifact.save()
def check_and_attempt_upload(
    artifact: Artifact,
    using_key: str | None = None,
    access_token: str | None = None,
    print_progress: bool = True,
    **kwargs,
) -> Exception | None:
    # kwargs are propagated to .upload_from in the end
    # if Artifact object is either newly instantiated or replace() was called on
    # a local env it will have a _local_filepath and needs to be uploaded
    if hasattr(artifact, "_local_filepath"):
        try:
            storage_path, cache_path = upload_artifact(
                artifact,
                using_key,
                access_token=access_token,
                print_progress=print_progress,
                **kwargs,
            )
        except Exception as exception:
            logger.warning(f"could not upload artifact: {artifact}")
            # clear dangling storages if we were actually uploading or saving
            if hasattr(artifact, "_to_store") and artifact._to_store:
                artifact._clear_storagekey = auto_storage_key_from_artifact(artifact)
            return exception
        # copies (if on-disk) or moves the temporary file (if in-memory) to the cache
        if os.getenv("LAMINDB_MULTI_INSTANCE") is None:
            # this happens only after the actual upload was performed
            # we avoid failing here in case any problems happen in copy_or_move_to_cache
            # because the cache copying or cleanup is not absolutely necessary
            try:
                copy_or_move_to_cache(artifact, storage_path, cache_path)
            except Exception as e:
                logger.warning(f"A problem with cache on saving: {e}")
        # after successful upload, we should remove the attribute so that another call
        # call to save won't upload again, the user should call replace() then
        del artifact._local_filepath
    # returning None means proceed (either success or no action needed)
    return None


def copy_or_move_to_cache(
    artifact: Artifact, storage_path: UPath, cache_path: UPath | None
):
    local_path = artifact._local_filepath

    # in-memory cases
    if local_path is None or not local_path.exists():
        return None

    local_path = local_path.resolve()
    is_dir = local_path.is_dir()
    cache_dir = settings.cache_dir

    # just delete from the cache dir if storage_path is local
    if cache_path is None:
        if (
            local_path.as_posix() != storage_path.as_posix()
            and cache_dir in local_path.parents
        ):
            if is_dir:
                shutil.rmtree(local_path)
            else:
                local_path.unlink()
        return None
    # non-local storage_path further
    if local_path != cache_path:
        if cache_path.exists():
            logger.warning(f"replacing the existing cache path {cache_path.as_posix()}")
            if cache_path.is_dir():
                shutil.rmtree(cache_path)
            else:
                cache_path.unlink()
        else:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
        if cache_dir in local_path.parents:
            local_path.replace(cache_path)
        else:
            if is_dir:
                shutil.copytree(local_path, cache_path)
            else:
                shutil.copy(local_path, cache_path)
    # make sure that the cached version is older than the cloud one
    mts = datetime.now().timestamp() + 1.0
    if is_dir:
        files = (file for file in cache_path.rglob("*") if file.is_file())
        for file in files:
            os.utime(file, times=(mts, mts))
    else:
        os.utime(cache_path, times=(mts, mts))


# This is also used within Artifact.save()
def check_and_attempt_clearing(
    artifact: Artifact,
    raise_file_not_found_error: bool = True,
    using_key: str | None = None,
) -> Exception | None:
    # this is a clean-up operation after replace() was called
    # or if there was an exception during upload
    if hasattr(artifact, "_clear_storagekey"):
        try:
            if artifact._clear_storagekey is not None:
                delete_msg = delete_storage_using_key(
                    artifact,
                    artifact._clear_storagekey,
                    raise_file_not_found_error=raise_file_not_found_error,
                    using_key=using_key,
                )
                if delete_msg != "did-not-delete":
                    logger.success(
                        f"deleted stale object at storage key {artifact._clear_storagekey}"
                    )
                artifact._clear_storagekey = None
        except Exception as exception:
            return exception
    # returning None means proceed (either success or no action needed)
    return None


def store_artifacts(
    artifacts: Iterable[Artifact], using_key: str | None = None
) -> None:
    """Upload artifacts in a list of database-committed artifacts to storage.

    If any upload fails, subsequent artifacts are cleaned up from the DB.
    """
    exception: Exception | None = None
    # because uploads might fail, we need to maintain a new list
    # of the succeeded uploads
    stored_artifacts = []

    # upload new local artifacts
    for artifact in artifacts:
        # failure here sets ._clear_storagekey
        # for cleanup below
        exception = check_and_attempt_upload(artifact, using_key)
        if exception is not None:
            break
        stored_artifacts += [artifact]
        # if check_and_attempt_upload was successfull
        # then this can have only ._clear_storagekey from .replace
        exception = check_and_attempt_clearing(
            artifact, raise_file_not_found_error=True, using_key=using_key
        )
        if exception is not None:
            logger.warning(f"clean up of {artifact._clear_storagekey} failed")
            break

    if exception is not None:
        # clean up metadata for artifacts not uploaded to storage
        with transaction.atomic():
            for artifact in artifacts:
                if artifact not in stored_artifacts:
                    artifact._delete_skip_storage()
                    # clean up storage after failure in check_and_attempt_upload
                    exception_clear = check_and_attempt_clearing(
                        artifact, raise_file_not_found_error=False, using_key=using_key
                    )
                    if exception_clear is not None:
                        logger.warning(
                            f"clean up of {artifact._clear_storagekey} after the upload error failed"
                        )
        error_message = prepare_error_message(artifacts, stored_artifacts, exception)
        # this is bad because we're losing the original traceback
        # needs to be refactored - also, the orginal error should be raised
        raise RuntimeError(error_message)
    return None


def prepare_error_message(records, stored_artifacts, exception) -> str:
    if len(stored_artifacts) == 0:
        error_message = (
            "No entries were uploaded or committed"
            " to the database. See error message:\n\n"
        )
    else:
        error_message = (
            "The following entries have been"
            " successfully uploaded and committed to the database:\n"
        )
        for record in stored_artifacts:
            error_message += (
                f"- {', '.join(record.__repr__().split(', ')[:3]) + ', ...)'}\n"
            )
        error_message += "\nSee error message:\n\n"
    error_message += f"{str(exception)}\n\n{traceback.format_exc()}"
    return error_message


def upload_artifact(
    artifact,
    using_key: str | None = None,
    access_token: str | None = None,
    print_progress: bool = True,
    **kwargs,
) -> tuple[UPath, UPath | None]:
    """Store and add file and its linked entries."""
    # kwargs are propagated to .upload_from in the end
    # can't currently use  filepath_from_artifact here because it resolves to ._local_filepath
    storage_key = auto_storage_key_from_artifact(artifact)
    storage_path, storage_settings = attempt_accessing_path(
        artifact, storage_key, using_key=using_key, access_token=access_token
    )
    if hasattr(artifact, "_to_store") and artifact._to_store:
        logger.save(f"storing artifact '{artifact.uid}' at '{storage_path}'")
        store_file_or_folder(
            artifact._local_filepath,
            storage_path,
            print_progress=print_progress,
            **kwargs,
        )

    if isinstance(storage_path, LocalPathClasses):
        cache_path = None
    else:
        cache_key = _cache_key_from_artifact_storage(artifact, storage_settings)
        cache_path = storage_settings.cloud_to_local_no_update(
            storage_path, cache_key=cache_key
        )

    return storage_path, cache_path
