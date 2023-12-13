import os
import shutil
import traceback
from collections import defaultdict
from datetime import datetime
from functools import partial
from typing import Iterable, List, Optional, Tuple, Union, overload

import lamindb_setup
from django.db import transaction
from django.utils.functional import partition
from lamin_utils import logger
from lamindb_setup.dev.upath import print_hook
from lnschema_core.models import Artifact, Registry

from lamindb.dev.storage.file import (
    auto_storage_key_from_artifact,
    delete_storage_using_key,
    store_artifact,
)

try:
    from lamindb.dev.storage._zarr import write_adata_zarr
except ImportError:

    def write_adata_zarr(filepath):  # type: ignore
        raise ImportError("Please install zarr: pip install zarr")


def save(
    records: Iterable[Registry], ignore_conflicts: Optional[bool] = False, **kwargs
) -> None:
    """Bulk save to registries & storage.

    Note:

        This is a much faster than saving records using ``record.save()``.

    Warning:

        Bulk saving neither automatically creates related records nor updates
        existing records! Use ``record.save()`` for these use cases.

    Args:
        records: Multiple :class:`~lamindb.dev.Registry` objects.
        ignore_conflicts: If ``True``, do not error if some records violate a
           unique or another constraint. However, it won't inplace update the id
           fields of records. If you need records with ids, you need to query
           them from the database.
        **kwargs: Get kwargs related to parents.

    Examples:

        Save a collection of records in one transaction, which is much faster
        than writing a loop over ``projects.save()``:

        >>> labels = [ln.ULabel(f"Label {i}") for i in range(10)]
        >>> ln.save(projects)

        For a single record, use ``record.save()``:

        >>> transform = ln.Transform(name="My pipeline")
        >>> transform.save()

        Update a single existing record:

        >>> transform = ln.filter(ln.Transform, uid="0Cb86EZj").one()
        >>> transform.name = "New name"
        >>> transform.save()

    """
    if isinstance(records, Registry):
        raise ValueError("Please use record.save() if saving a single record.")

    # previously, this was all set based,
    # but models without primary keys aren't hashable
    # we distinguish between artifacts and non-artifacts
    # for artifacts, we want to bulk-upload rather than upload one-by-one
    non_artifacts, artifacts = partition(lambda r: isinstance(r, Artifact), records)
    if non_artifacts:
        # first save all records that do not yet have a primary key without
        # recursing parents
        _, non_artifacts_without_pk = partition(lambda r: r.pk is None, non_artifacts)
        bulk_create(non_artifacts_without_pk, ignore_conflicts=ignore_conflicts)
        non_artifacts_with_parents = [
            r for r in non_artifacts_without_pk if hasattr(r, "_parents")
        ]

        if len(non_artifacts_with_parents) > 0 and kwargs.get("parents") is not False:
            # this can only happen within lnschema_bionty right now!!
            # we might extend to core lamindb later
            import lnschema_bionty as lb

            if kwargs.get("parents") or (
                kwargs.get("parents") is None and lb.settings.auto_save_parents
            ):
                mute = False if kwargs.get("mute") is None else kwargs.get("mute")
                if not mute:
                    # save the record with parents one by one
                    logger.warning(
                        "now recursing through parents: "
                        "this only happens once, but is much slower than bulk saving"
                    )
                    logger.hint(
                        "you can switch this off via: lb.settings.auto_save_parents ="
                        " False"
                    )
                for record in non_artifacts_with_parents:
                    record._save_ontology_parents(mute=True)

    if artifacts:
        with transaction.atomic():
            for record in artifacts:
                record._save_skip_storage()
        store_artifacts(artifacts)

    # this function returns None as potentially 10k records might be saved
    # refreshing all of them from the DB would mean a severe performance penalty
    # 2nd reason: consistency with Django Model.save(), which also returns None
    return None


def bulk_create(records: Iterable[Registry], ignore_conflicts: Optional[bool] = False):
    records_by_orm = defaultdict(list)
    for record in records:
        records_by_orm[record.__class__].append(record)
    for orm, records in records_by_orm.items():
        orm.objects.bulk_create(records, ignore_conflicts=ignore_conflicts)


# This is also used within Artifact.save()
def check_and_attempt_upload(artifact: Artifact) -> Optional[Exception]:
    # if Artifact object is either newly instantiated or replace() was called on
    # a local env it will have a _local_filepath and needs to be uploaded
    if hasattr(artifact, "_local_filepath"):
        try:
            upload_artifact(artifact)
        except Exception as exception:
            logger.warning(f"could not upload artifact: {artifact}")
            return exception
        # copies (if on-disk) or moves the temporary file (if in-memory) to the cache
        copy_or_move_to_cache(artifact)
        # after successful upload, we should remove the attribute so that another call
        # call to save won't upload again, the user should call replace() then
        del artifact._local_filepath
    # returning None means proceed (either success or no action needed)
    return None


def copy_or_move_to_cache(artifact: Artifact):
    local_path = artifact._local_filepath

    # in-memory zarr or on-disk zarr
    if local_path is None or not local_path.is_file():
        return None

    local_path = local_path.resolve()
    cache_dir = lamindb_setup.settings.storage.cache_dir

    # local instance, just delete the cached file
    if not lamindb_setup.settings.storage.is_cloud:
        if cache_dir in local_path.parents:
            local_path.unlink()
        return None

    # maybe create something like storage.key_to_local(key) later to simplfy
    storage_key = auto_storage_key_from_artifact(artifact)
    storage_path = lamindb_setup.settings.storage.key_to_filepath(storage_key)
    cache_path = lamindb_setup.settings.storage.cloud_to_local_no_update(storage_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_dir in local_path.parents:
        local_path.replace(cache_path)
    else:
        shutil.copy(local_path, cache_path)
    # make sure that the cached version is older than the cloud one
    mts = datetime.now().timestamp() + 1.0
    os.utime(cache_path, times=(mts, mts))


# This is also used within Artifact.save()
def check_and_attempt_clearing(artifact: Artifact) -> Optional[Exception]:
    # this is a clean-up operation after replace() was called
    # this will only evaluate to True if replace() was called
    if hasattr(artifact, "_clear_storagekey"):
        try:
            if artifact._clear_storagekey is not None:
                delete_storage_using_key(artifact, artifact._clear_storagekey)
                logger.success(
                    f"deleted stale object at storage key {artifact._clear_storagekey}"
                )
                artifact._clear_storagekey = None
        except Exception as exception:
            return exception
    # returning None means proceed (either success or no action needed)
    return None


def store_artifacts(artifacts: Iterable[Artifact]) -> None:
    """Upload artifacts in a list of database-committed artifacts to storage.

    If any upload fails, subsequent artifacts are cleaned up from the DB.
    """
    exception: Optional[Exception] = None
    # because uploads might fail, we need to maintain a new list
    # of the succeeded uploads
    stored_artifacts = []

    # upload new local artifacts
    for artifact in artifacts:
        exception = check_and_attempt_upload(artifact)
        if exception is not None:
            break
        stored_artifacts += [artifact]
        exception = check_and_attempt_clearing(artifact)
        if exception is not None:
            logger.warning(f"clean up of {artifact._clear_storagekey} failed")
            break

    if exception is not None:
        # clean up metadata for artifacts not uploaded to storage
        with transaction.atomic():
            for artifact in artifacts:
                if artifact not in stored_artifacts:
                    artifact._delete_skip_storage()
        error_message = prepare_error_message(artifacts, stored_artifacts, exception)
        raise RuntimeError(error_message)
    return None


def prepare_error_message(records, stored_artifacts, exception) -> str:
    if len(records) == 1 or len(stored_artifacts) == 0:
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


def upload_artifact(artifact) -> None:
    """Store and add file and its linked entries."""
    # do NOT hand-craft the storage key!
    artifact_storage_key = auto_storage_key_from_artifact(artifact)
    storage_path = lamindb_setup.settings.instance.storage.key_to_filepath(
        artifact_storage_key
    )
    msg = f"storing artifact '{artifact.uid}' at '{storage_path}'"
    if (
        artifact.suffix in {".zarr", ".zrad"}
        and hasattr(artifact, "_memory_rep")
        and artifact._memory_rep is not None
    ):
        logger.save(msg)
        print_progress = partial(
            print_hook, filepath=artifact_storage_key, action="uploading"
        )
        write_adata_zarr(artifact._memory_rep, storage_path, callback=print_progress)
    elif hasattr(artifact, "_to_store") and artifact._to_store:
        logger.save(msg)
        store_artifact(artifact._local_filepath, artifact_storage_key)
