# ruff: noqa: TC004
from __future__ import annotations

import shutil
import warnings
from collections import defaultdict
from pathlib import Path, PurePath, PurePosixPath
from typing import TYPE_CHECKING, Any, Literal, Union, overload

import fsspec
import lamindb_setup as ln_setup
import pandas as pd
from anndata import AnnData
from django.db import ProgrammingError, models
from django.db.models import CASCADE, PROTECT, Q
from django.db.models.functions import Length
from lamin_utils import colors, logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core._hub_core import select_storage_or_parent
from lamindb_setup.core.hashing import HASH_LENGTH, hash_dir, hash_file
from lamindb_setup.core.upath import (
    create_path,
    extract_suffix_from_path,
    get_stat_dir_cloud,
    get_stat_file_cloud,
)
from lamindb_setup.types import UPathStr

from lamindb.base import deprecated
from lamindb.base.fields import (
    BigIntegerField,
    BooleanField,
    CharField,
    ForeignKey,
    TextField,
)
from lamindb.errors import FieldValidationError, NoWriteAccess, UnknownStorageLocation
from lamindb.models.query_set import QuerySet, SQLRecordList

from ..base.users import current_user_id
from ..core._settings import is_read_only_connection, settings
from ..core.loaders import load_to_memory
from ..core.storage import (
    LocalPathClasses,
    UPath,
    delete_storage,
    infer_suffix,
    write_to_disk,
)
from ..core.storage._anndata_accessor import _anndata_n_observations
from ..core.storage._backed_access import (
    _track_writes_factory,
    backed_access,
)
from ..core.storage._polars_lazy_df import POLARS_SUFFIXES
from ..core.storage._pyarrow_dataset import PYARROW_SUFFIXES
from ..core.storage._tiledbsoma import _soma_n_observations
from ..core.storage.paths import (
    AUTO_KEY_PREFIX,
    auto_storage_key_from_artifact,
    auto_storage_key_from_artifact_uid,
    check_path_is_child_of_root,
    filepath_cache_key_from_artifact,
    filepath_from_artifact,
)
from ..errors import InvalidArgument, NoStorageLocationForSpace, ValidationError
from ..models._is_versioned import (
    create_uid,
)
from ._feature_manager import (
    FeatureManager,
    get_label_links,
)
from ._is_versioned import IsVersioned
from ._relations import (
    dict_module_name_to_model_name,
    dict_related_model_to_related_name,
)
from .feature import Feature, FeatureValue
from .has_parents import view_lineage
from .run import Run, TracksRun, TracksUpdates, User
from .save import check_and_attempt_clearing, check_and_attempt_upload
from .schema import Schema
from .sqlrecord import (
    BaseSQLRecord,
    IsLink,
    SQLRecord,
    _get_record_kwargs,
)
from .storage import Storage
from .ulabel import ULabel

WARNING_RUN_TRANSFORM = "no run & transform got linked, call `ln.track()` & re-run"

WARNING_NO_INPUT = "run input wasn't tracked, call `ln.track()` and re-run"

try:
    from ..core.storage._zarr import identify_zarr_type
except ImportError:

    def identify_zarr_type(storepath):  # type: ignore
        raise ImportError("Please install zarr: pip install 'lamindb[zarr]'")


if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from mudata import MuData  # noqa: TC004
    from polars import LazyFrame as PolarsLazyFrame
    from pyarrow.dataset import Dataset as PyArrowDataset
    from spatialdata import SpatialData  # noqa: TC004
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from tiledbsoma import Measurement as SOMAMeasurement

    from lamindb.base.types import StrField
    from lamindb.core.storage._backed_access import (
        AnnDataAccessor,
        BackedAccessor,
        SpatialDataAccessor,
    )
    from lamindb.core.types import ScverseDataStructures

    from ..base.types import (
        ArtifactKind,
    )
    from ._label_manager import LabelManager
    from .block import ArtifactBlock
    from .collection import Collection
    from .project import Project, Reference
    from .record import Record
    from .sqlrecord import Branch, Space
    from .transform import Transform


INCONSISTENT_STATE_MSG = (
    "Trying to read a folder artifact from an outdated version, "
    "this can result in an incosistent state.\n"
    "Read from the latest version: artifact.versions.get(is_latest=True)"
)


def process_pathlike(
    filepath: UPath,
    storage: Storage,
    using_key: str | None,
    skip_existence_check: bool = False,
) -> tuple[Storage, bool]:
    """Determines the appropriate storage for a given path and whether to use an existing storage key."""
    if not skip_existence_check:
        try:  # check if file exists
            if not filepath.exists():
                raise FileNotFoundError(filepath)
        except PermissionError:
            pass
    if check_path_is_child_of_root(filepath, storage.root):
        use_existing_storage_key = True
        return storage, use_existing_storage_key
    else:
        # check whether the path is part of one of the existing
        # already-registered storage locations
        result = None
        # within the hub, we don't want to perform check_path_in_existing_storage
        if using_key is None:
            result = check_path_in_existing_storage(
                filepath, check_hub_register_storage=setup_settings.instance.is_on_hub
            )
        if isinstance(result, Storage):
            use_existing_storage_key = True
            return result, use_existing_storage_key
        else:
            # if the path is in the cloud, we have a good candidate
            # for the storage root: the bucket
            if not isinstance(filepath, LocalPathClasses):
                # for a cloud path, new_root is always the bucket name
                if filepath.protocol == "hf":
                    hf_path = filepath.fs.resolve_path(filepath.as_posix())
                    hf_path.path_in_repo = ""
                    new_root = "hf://" + hf_path.unresolve().rstrip("/")
                else:
                    if filepath.protocol == "s3":
                        # check that endpoint_url didn't propagate here
                        # as a part of the path string
                        assert "?" not in filepath.path  # noqa: S101
                    new_root = list(filepath.parents)[-1].as_posix().rstrip("/")
                # Re the Parallel execution of the logic below:
                # One of the threads (or processes) would start to write the hub record and then the test file.
                # The other ones would retrieve the hub record and the test file.
                # All of them would come out of the exercise with storage_record.instance_uid == setup_settings.instance.uid
                # and all of them would raise UnkownStorageLocation.
                # Then one of these threads will trigger storage_record.delete() but also this is idempotent;
                # this means they all throw the same error and deletion of the inexistent stuff (hub record, marker file)
                # would just silently fail.
                # Edge case: A user legitimately creates a storage location and another user runs this here at the exact same time.
                # There is no way to decide then which is the legitimate creation.
                storage_record = Storage(root=new_root).save()
                if storage_record.instance_uid == setup_settings.instance.uid:
                    # we don't want to inadvertently create managed storage locations
                    # hence, we revert the creation and throw an error
                    storage_record.delete()
                    raise UnknownStorageLocation(
                        f"Path {filepath} is not contained in any known storage location:\n{Storage.to_dataframe()[['uid', 'root', 'type']]}\n\n"
                        f"Create a managed storage location that contains the path, e.g., by calling: ln.Storage(root='{new_root}').save()"
                    )
                use_existing_storage_key = True
                return storage_record, use_existing_storage_key
            # if the filepath is local
            else:
                use_existing_storage_key = False
                # if the default storage is local we'll throw an error if the user
                # doesn't provide a key
                if storage.type == "local":
                    return storage, use_existing_storage_key
                # if the default storage is in the cloud (the file is going to
                # be uploaded upon saving it), we treat the filepath as a cache
                else:
                    return storage, use_existing_storage_key


def process_data(
    provisional_uid: str,
    data: UPathStr | pd.DataFrame | AnnData,
    format: str | None,
    key: str | None,
    storage: Storage,
    using_key: str | None,
    skip_existence_check: bool = False,
    is_replace: bool = False,
    to_disk_kwargs: dict[str, Any] | None = None,
) -> tuple[Any, Path | UPath, str, Storage, bool]:
    """Serialize a data object that's provided as file or in memory.

    if not overwritten, data gets stored in default storage
    """
    if key is not None:
        key_suffix = extract_suffix_from_path(PurePosixPath(key), arg_name="key")
        # use suffix as the (adata) format if the format is not provided
        if isinstance(data, AnnData) and format is None and len(key_suffix) > 0:
            format = key_suffix[1:]
    else:
        key_suffix = None

    if isinstance(data, (str, Path, UPath)):
        access_token = (
            storage._access_token if hasattr(storage, "_access_token") else None
        )
        path = create_path(data, access_token=access_token)
        # we don't resolve http links because they can resolve into a different domain
        # for example into a temporary url
        if path.protocol not in {"http", "https"}:
            path = path.resolve()

        storage, use_existing_storage_key = process_pathlike(
            path,
            storage=storage,
            using_key=using_key,
            skip_existence_check=skip_existence_check,
        )
        suffix = extract_suffix_from_path(path)
        memory_rep = None
    elif (
        isinstance(data, pd.DataFrame)
        or isinstance(data, AnnData)
        or data_is_scversedatastructure(data, "MuData")
        or data_is_scversedatastructure(data, "SpatialData")
    ):
        storage = storage
        memory_rep = data
        suffix = infer_suffix(data, format)
    else:
        raise NotImplementedError(
            f"Do not know how to create an Artifact from {data}, pass a path instead."
        )

    # Check for suffix consistency
    if key_suffix is not None and key_suffix != suffix and not is_replace:
        # consciously omitting a trailing period
        if isinstance(data, (str, Path, UPath)):  # UPathStr, spelled out
            message = f"The passed path's suffix '{suffix}' must match the passed key's suffix '{key_suffix}'."
        else:
            message = f"The passed key's suffix '{key_suffix}' must match the passed path's suffix '{suffix}'."
        raise InvalidArgument(message)

    # in case we have an in-memory representation, we need to write it to disk
    if memory_rep is not None:
        path = settings.cache_dir / f"{provisional_uid}{suffix}"
        logger.important("writing the in-memory object into cache")
        if to_disk_kwargs is None:
            to_disk_kwargs = {}
        write_to_disk(data, path, **to_disk_kwargs)
        use_existing_storage_key = False

    return memory_rep, path, suffix, storage, use_existing_storage_key


def get_stat_or_artifact(
    path: UPath,
    storage: Record,
    key: str | None = None,
    check_hash: bool = True,
    is_replace: bool = False,
    instance: str | None = None,
    skip_hash_lookup: bool = False,
) -> Union[tuple[int, str | None, str | None, int | None, Artifact | None], Artifact]:
    """Retrieves file statistics or an existing artifact based on the path, hash, and key."""
    n_files = None
    if settings.creation.artifact_skip_size_hash:
        return None, None, None, n_files, None
    stat = path.stat()  # one network request
    if not isinstance(path, LocalPathClasses):
        size, hash, hash_type = None, None, None
        if stat is not None:
            # convert UPathStatResult to fsspec info dict
            stat = stat.as_info()
            if (store_type := stat["type"]) == "file":
                size, hash, hash_type = get_stat_file_cloud(stat)
            elif store_type == "directory":
                size, hash, hash_type, n_files = get_stat_dir_cloud(path)
        if hash is None:
            logger.warning(f"did not add hash for {path}")
            return size, hash, hash_type, n_files, None
    else:
        if path.is_dir():
            size, hash, hash_type, n_files = hash_dir(path)
        else:
            hash, hash_type = hash_file(path)
            size = stat.st_size
    if not check_hash:
        return size, hash, hash_type, n_files, None
    previous_artifact_version = None
    artifacts_qs = Artifact.objects.using(instance)
    if skip_hash_lookup:
        artifact_with_same_hash_exists = False
        if key is not None and not is_replace:
            # only search for a previous version of the artifact
            # ignoring hash
            lookup_result = artifacts_qs.filter(
                ~Q(branch_id=-1),
                key=key,
                storage=storage,
            ).order_by("-created_at")
        else:
            lookup_result = []
    else:
        # this purposefully leaves out the storage location and key that we have
        # in the hard database unique constraints
        # so that the user is able to find artifacts with the same hash across
        # storage locations and keys
        # if this is not desired, set skip_hash_lookup=True
        if key is None or is_replace:
            lookup_result = artifacts_qs.filter(~Q(branch_id=-1), hash=hash)
            artifact_with_same_hash_exists = len(lookup_result) > 0
        else:
            # the following query achieves one more thing beyond hash lookup
            # it allows us to find a previous version of the artifact based on
            # matching key & storage even if the hash is different
            # we do this here so that we don't have to do an additional query later
            # see the `previous_artifact_version` variable below
            lookup_result = artifacts_qs.filter(
                ~Q(branch_id=-1),
                Q(hash=hash) | Q(key=key, storage=storage),
            ).order_by("-created_at")
            artifact_with_same_hash_exists = lookup_result.filter(hash=hash).count() > 0
    if key is not None and not is_replace:
        if not artifact_with_same_hash_exists and len(lookup_result) > 0:
            logger.important(
                f"creating new artifact version for key '{key}' in storage '{storage.root}'"
            )
            previous_artifact_version = lookup_result[0]
    if artifact_with_same_hash_exists:
        message = "returning artifact with same hash"
        logger.important(
            f"{message}: {lookup_result[0]}; to track this artifact as an input, use: ln.Artifact.get()"
        )
        return lookup_result[0]
    else:
        return size, hash, hash_type, n_files, previous_artifact_version


def check_path_in_existing_storage(
    path: Path | UPath,
    check_hub_register_storage: bool = False,
    using_key: str | None = None,
) -> Storage | None:
    for storage in Storage.objects.using(using_key).order_by(Length("root").desc()):
        # if path is part of storage, return it
        if check_path_is_child_of_root(path, root=storage.root):
            return storage
    # we don't see parents registered in the db, so checking the hub
    # just check for 2 writable cloud protocols, maybe change in the future
    if check_hub_register_storage and getattr(path, "protocol", None) in {"s3", "gs"}:
        result = select_storage_or_parent(path.as_posix())
        if result is not None:
            return Storage(**result, _skip_preparation=True).save()
    return None


def get_relative_path_to_directory(
    path: PurePath | Path | UPath, directory: PurePath | Path | UPath
) -> PurePath | Path:
    if isinstance(directory, UPath) and not isinstance(directory, LocalPathClasses):
        # UPath.relative_to() is not behaving as it should (2023-04-07)
        # need to lstrip otherwise inconsistent behavior across trailing slashes
        # see test_artifact.py: test_get_relative_path_to_directory
        relpath = PurePath(
            path.as_posix().replace(directory.as_posix(), "").lstrip("/")
        )
    elif isinstance(directory, Path):
        relpath = path.resolve().relative_to(directory.resolve())  # type: ignore
    elif isinstance(directory, PurePath):
        relpath = path.relative_to(directory)
    else:
        raise TypeError("Directory not of type Path or UPath")
    return relpath


def get_artifact_kwargs_from_data(
    *,
    data: Path | UPath | str | pd.DataFrame | ScverseDataStructures,
    key: str | None,
    run: Run | None,
    format: str | None,
    provisional_uid: str,
    version: str | None,
    storage: Storage,
    using_key: str | None = None,
    is_replace: bool = False,
    skip_check_exists: bool = False,
    overwrite_versions: bool | None = None,
    skip_hash_lookup: bool = False,
    to_disk_kwargs: dict[str, Any] | None = None,
):
    memory_rep, path, suffix, storage, use_existing_storage_key = process_data(
        provisional_uid,
        data,
        format,
        key,
        storage,
        using_key,
        skip_check_exists,
        is_replace=is_replace,
        to_disk_kwargs=to_disk_kwargs,
    )

    check_path_in_storage = False
    real_key = None
    if use_existing_storage_key:
        inferred_key = get_relative_path_to_directory(
            path=path, directory=UPath(storage.root)
        ).as_posix()
        if key is None:
            key = inferred_key
        elif key != inferred_key:
            real_key = inferred_key
        check_path_in_storage = True
    else:
        storage = storage

    stat_or_artifact = get_stat_or_artifact(
        path=path,
        storage=storage,
        key=key,
        instance=using_key,
        is_replace=is_replace,
        skip_hash_lookup=skip_hash_lookup,
    )
    if isinstance(stat_or_artifact, Artifact):
        existing_artifact = stat_or_artifact
        return existing_artifact, key
    else:
        size, hash, hash_type, n_files, revises = stat_or_artifact

    if revises is not None:  # update provisional_uid
        provisional_uid, revises = create_uid(revises=revises, version=version)
        if settings.cache_dir in path.parents:
            path = path.rename(path.with_name(f"{provisional_uid}{suffix}"))

    log_storage_hint(
        check_path_in_storage=check_path_in_storage,
        storage=storage,
        key=key,
        uid=provisional_uid,
        suffix=suffix,
        is_dir=n_files is not None,
    )

    if overwrite_versions is None:
        overwrite_versions = n_files is not None

    if check_path_in_storage:
        # we use an actual storage key if key is not provided explicitly
        key_is_virtual = real_key is not None
    else:
        # do we use a virtual or an actual storage key?
        key_is_virtual = settings.creation._artifact_use_virtual_keys

    kwargs = {
        "uid": provisional_uid,
        "suffix": suffix,
        "hash": hash,
        "_hash_type": hash_type,
        "key": key,
        "size": size,
        "storage_id": storage.id,
        # passing both the id and the object
        # to make them both available immediately
        # after object creation
        "n_files": n_files,
        "_overwrite_versions": overwrite_versions,  # True for folder, False for file
        "n_observations": None,  # to implement
        "run_id": run.id if run is not None else None,
        "run": run,
        "_key_is_virtual": key_is_virtual,
        "revises": revises,
        "_real_key": real_key,
    }
    if not isinstance(path, LocalPathClasses):
        local_filepath = None
        cloud_filepath = path
    else:
        local_filepath = path
        cloud_filepath = None
    privates = {
        "local_filepath": local_filepath,
        "cloud_filepath": cloud_filepath,
        "memory_rep": memory_rep,
        "check_path_in_storage": check_path_in_storage,
    }
    return kwargs, privates


def log_storage_hint(
    *,
    check_path_in_storage: bool,
    storage: Storage | None,
    key: str | None,
    uid: str,
    suffix: str,
    is_dir: bool,
) -> None:
    hint = ""
    if check_path_in_storage:
        display_root = storage.root  # type: ignore
        # check whether path is local
        if fsspec.utils.get_protocol(storage.root) == "file":  # type: ignore
            # if it's a local path, check whether it's in the current working directory
            root_path = Path(storage.root)  # type: ignore
            if check_path_is_child_of_root(root_path, Path.cwd()):
                # only display the relative path, not the fully resolved path
                display_root = root_path.relative_to(Path.cwd())  # type: ignore
        hint += f"path in storage '{display_root}'"  # type: ignore
    else:
        hint += "path content will be copied to default storage upon `save()`"
    if key is None:
        storage_key = auto_storage_key_from_artifact_uid(uid, suffix, is_dir)
        hint += f" with key `None` ('{storage_key}')"
    else:
        hint += f" with key '{key}'"
    logger.hint(hint)


def data_is_scversedatastructure(
    data: ScverseDataStructures | UPathStr,
    structure_type: Literal["AnnData", "MuData", "SpatialData"] | None = None,
    cloud_warning: bool = True,
) -> bool:
    """Determine whether a specific in-memory object or a UPathstr is any or a specific scverse data structure."""
    file_suffix = None
    if structure_type == "AnnData":
        file_suffix = ".h5ad"
    elif structure_type == "MuData":
        file_suffix = ".h5mu"
    # SpatialData does not have a unique suffix but `.zarr`

    # AnnData allows both AnnDataAccessor and AnnData
    class_name = data.__class__.__name__
    if structure_type is None:
        return any(
            class_name
            in (["AnnData", "AnnDataAccessor"] if cl_name == "AnnData" else [cl_name])
            for cl_name in ["AnnData", "MuData", "SpatialData"]
        )
    elif class_name in (
        ["AnnData", "AnnDataAccessor"]
        if structure_type == "AnnData"
        else [structure_type]
    ):
        return True

    data_type = structure_type.lower()
    if isinstance(data, (str, Path, UPath)):
        data_path = UPath(data)

        if file_suffix in data_path.suffixes:
            return True

        if data_path.suffix == ".zarr":
            type_suffix = f".{data_type}"
            if type_suffix in data_path.suffixes:
                return True

            # check only for local, expensive for cloud
            if fsspec.utils.get_protocol(data_path.as_posix()) == "file":
                return (
                    identify_zarr_type(
                        data_path if structure_type == "AnnData" else data,
                        check=True if structure_type == "AnnData" else False,
                    )
                    == data_type
                )
            elif cloud_warning:
                logger.warning(
                    f"we do not check whether cloud zarr is {structure_type}"
                )
                return False

    return False


def data_is_soma_experiment(data: SOMAExperiment | UPathStr) -> bool:
    # We are not importing tiledbsoma here to keep loaded modules minimal
    if hasattr(data, "__class__") and data.__class__.__name__ == "Experiment":
        return True
    if isinstance(data, (str, Path)):
        return UPath(data).suffix == ".tiledbsoma"
    return False


def _check_otype_artifact(
    data: UPathStr | pd.DataFrame | ScverseDataStructures,
    otype: str | None = None,
    cloud_warning: bool = True,
) -> str:
    if otype is None:
        if isinstance(data, pd.DataFrame):
            logger.warning("data is a DataFrame, please use .from_dataframe()")
            otype = "DataFrame"
            return otype

        data_is_path = isinstance(data, (str, Path))
        if data_is_scversedatastructure(data, "AnnData", cloud_warning):
            if not data_is_path:
                logger.warning("data is an AnnData, please use .from_anndata()")
            otype = "AnnData"
        elif data_is_scversedatastructure(data, "MuData", cloud_warning):
            if not data_is_path:
                logger.warning("data is a MuData, please use .from_mudata()")
            otype = "MuData"
        elif data_is_scversedatastructure(data, "SpatialData", cloud_warning):
            if not data_is_path:
                logger.warning("data is a SpatialData, please use .from_spatialdata()")
            otype = "SpatialData"
        elif not data_is_path:  # UPath is a subclass of Path
            raise TypeError("data has to be a string, Path, UPath")
    return otype


def populate_subsequent_run(record: Artifact | Collection, run: Run | None) -> None:
    if run is None:
        return
    if record.run is None:
        record.run = run
    elif record.run != run:
        record._subsequent_runs.add(run)
        record._subsequent_run_id = run.id


# also see current_run() in core._data
def get_run(run: Run | None) -> Run | None:
    from .._tracked import get_current_tracked_run
    from ..core._context import context

    if run is None:
        run = get_current_tracked_run()
        if run is None:
            run = context.run
        if run is None and not settings.creation.artifact_silence_missing_run_warning:
            if not is_read_only_connection():
                logger.warning(WARNING_RUN_TRANSFORM)
    # suppress run by passing False
    elif not run:
        run = None
    return run


def save_staged_feature_sets(self: Artifact) -> None:
    if hasattr(self, "_staged_feature_sets"):
        from lamindb.models._feature_manager import get_schema_by_slot_

        existing_staged_feature_sets = get_schema_by_slot_(self)
        saved_staged_feature_sets = {}
        for key, schema in self._staged_feature_sets.items():
            if isinstance(schema, Schema) and schema._state.adding:
                schema.save()
                saved_staged_feature_sets[key] = schema
            if key in existing_staged_feature_sets:
                # remove existing feature set on the same slot
                self.feature_sets.remove(existing_staged_feature_sets[key])
        if len(saved_staged_feature_sets) > 0:
            s = "s" if len(saved_staged_feature_sets) > 1 else ""
            display_schema_keys = ",".join(
                f"'{key}'" for key in saved_staged_feature_sets.keys()
            )
            logger.save(
                f"saved {len(saved_staged_feature_sets)} feature set{s} for slot{s}:"
                f" {display_schema_keys}"
            )


def save_schema_links(self: Artifact) -> None:
    from lamindb.models.save import bulk_create

    if hasattr(self, "_staged_feature_sets"):
        links = []
        for slot, schema in self._staged_feature_sets.items():
            kwargs = {
                "artifact_id": self.id,
                "schema_id": schema.id,
                "slot": slot,
            }
            links.append(Artifact.feature_sets.through(**kwargs))
        bulk_create(links, ignore_conflicts=True)


def validate_feature(feature: Feature, records: list[SQLRecord]) -> None:
    """Validate feature record, adjust feature.dtype based on labels records."""
    if not isinstance(feature, Feature):
        raise TypeError("feature has to be of type Feature")
    if feature._state.adding:
        registries = {record.__class__.__get_name_with_module__() for record in records}
        registries_str = "|".join(registries)
        msg = f"ln.Feature(name='{feature.name}', type='cat[{registries_str}]').save()"
        raise ValidationError(f"Feature not validated. If it looks correct: {msg}")


def get_labels(
    self,
    feature: Feature,
    mute: bool = False,
    flat_names: bool = False,
) -> QuerySet | dict[str, QuerySet] | list:
    """{}"""  # noqa: D415
    from .record import Record

    if not isinstance(feature, Feature):
        raise TypeError("feature has to be of type Feature")
    if feature.dtype is None or not feature.dtype.startswith("cat["):
        raise ValueError("feature does not have linked labels")
    registries_to_check = feature.dtype.replace("cat[", "").rstrip("]").split("|")
    if len(registries_to_check) > 1 and not mute:
        logger.warning("labels come from multiple registries!")
    # return an empty query set if self.id is still None
    if self.id is None:
        return QuerySet(self.__class__)
    qs_by_registry = {}
    for registry in registries_to_check:
        # currently need to distinguish between ULabel and non-ULabel, because
        # we only have the feature information for Label
        if registry in {"ULabel", "Record"}:
            links_to_labels = get_label_links(self, registry, feature)
            label_ids = [
                (link.ulabel_id if registry == "ULabel" else link.record_id)
                for link in links_to_labels
            ]
            model = ULabel if registry == "ULabel" else Record
            qs_by_registry[registry] = model.objects.using(self._state.db).filter(
                id__in=label_ids
            )
        elif registry in self.features._accessor_by_registry:
            qs_by_registry[registry] = getattr(
                self, self.features._accessor_by_registry[registry]
            ).all()
    if flat_names:
        # returns a flat list of names
        from .sqlrecord import get_name_field

        values = []
        for v in qs_by_registry.values():
            values += v.to_list(get_name_field(v))
        return values
    if len(registries_to_check) == 1 and registry in qs_by_registry:
        return qs_by_registry[registry]
    else:
        return qs_by_registry


def add_labels(
    self,
    records: SQLRecord | list[SQLRecord] | QuerySet | Iterable,
    feature: Feature | None = None,
    *,
    field: StrField | None = None,
    feature_ref_is_name: bool | None = None,
    label_ref_is_name: bool | None = None,
    from_curator: bool = False,
) -> None:
    """{}"""  # noqa: D415
    if self._state.adding:
        raise ValueError("Please save the artifact/collection before adding a label!")

    if isinstance(records, (QuerySet, QuerySet.__base__)):  # need to have both
        records = records.to_list()
    if isinstance(records, (str, SQLRecord)):
        records = [records]
    if not isinstance(records, list):  # avoids warning for pd Series
        records = list(records)
    # create records from values
    if len(records) == 0:
        return None
    if isinstance(records[0], str):  # type: ignore
        records_validated = []
        # feature is needed if we want to create records from values
        if feature is None:
            raise ValueError(
                "Please pass a feature, e.g., via: label = ln.ULabel(name='my_label',"
                " feature=ln.Feature(name='my_feature'))"
            )
        if feature.dtype.startswith("cat["):
            orm_dict = dict_module_name_to_model_name(Artifact)
            for reg in feature.dtype.replace("cat[", "").rstrip("]").split("|"):
                registry = orm_dict.get(reg)
                records_validated += registry.from_values(records, field=field)

        # feature doesn't have registries and therefore can't create records from values
        # ask users to pass records
        if len(records_validated) == 0:
            raise ValueError(
                "Please pass a record (a `SQLRecord` object), not a string, e.g., via:"
                " label"
                f" = ln.Record(name='{records[0]}')"  # type: ignore
            )
        records = records_validated

    for record in records:
        if record._state.adding:
            raise ValidationError(
                f"{record} not validated. If it looks correct: record.save()"
            )

    if feature is None:
        d = dict_related_model_to_related_name(self.__class__)
        # strategy: group records by registry to reduce number of transactions
        records_by_related_name: dict = {}
        for record in records:
            related_name = d.get(record.__class__.__get_name_with_module__())
            if related_name is None:
                raise ValueError(f"Can't add labels to {record.__class__} record!")
            if related_name not in records_by_related_name:
                records_by_related_name[related_name] = []
            records_by_related_name[related_name].append(record)
        for related_name, records in records_by_related_name.items():
            getattr(self, related_name).add(*records)
    else:
        validate_feature(feature, records)  # type:ignore
        records_by_registry = defaultdict(list)
        feature_sets = self.feature_sets.filter(itype="Feature")
        internal_features = set()  # type: ignore
        if len(feature_sets) > 0:
            for schema in feature_sets:
                internal_features = internal_features.union(
                    set(schema.members.values_list("name", flat=True))
                )  # type: ignore
        for record in records:
            records_by_registry[record.__class__.__get_name_with_module__()].append(
                record
            )
        for registry_name, records in records_by_registry.items():
            if not from_curator and feature.name in internal_features:
                raise ValidationError(
                    "Cannot manually annotate a feature measured *within* the dataset. Please use a Curator."
                )
            if registry_name not in feature.dtype:
                if not feature.dtype.startswith("cat"):
                    raise ValidationError(
                        f"Feature {feature.name} needs dtype='cat' for label annotation, currently has dtype='{feature.dtype}'"
                    )
                if feature.dtype == "cat":
                    feature.dtype = f"cat[{registry_name}]"  # type: ignore
                    feature.save()
                elif registry_name not in feature.dtype:
                    new_dtype = feature.dtype.rstrip("]") + f"|{registry_name}]"
                    raise ValidationError(
                        f"Label type {registry_name} is not valid for Feature(name='{feature.name}', dtype='{feature.dtype}'), consider updating to dtype='{new_dtype}'"
                    )

            if registry_name not in self.features._accessor_by_registry:
                logger.warning(f"skipping {registry_name}")
                continue
            if len(records) == 0:
                continue
            features_labels = {
                registry_name: [(feature, label_record) for label_record in records]
            }
            self.features._add_label_feature_links(
                features_labels,
                feature_ref_is_name=feature_ref_is_name,
                label_ref_is_name=label_ref_is_name,
            )


def delete_permanently(artifact: Artifact, storage: bool, using_key: str):
    # need to grab file path before deletion
    try:
        path, _ = filepath_from_artifact(artifact, using_key)
    except OSError:
        # we can still delete the record
        logger.warning("Could not get path")
        storage = False
    # only delete in storage if DB delete is successful
    # DB delete might error because of a foreign key constraint violated etc.
    if artifact._overwrite_versions and artifact.is_latest:
        logger.important(
            "deleting all versions of this artifact because they all share the same store"
        )
        for version in artifact.versions.all():  # includes artifact
            _delete_skip_storage(version)
    else:
        artifact._delete_skip_storage()
    # by default do not delete storage if deleting only a previous version
    # and the underlying store is mutable
    if artifact._overwrite_versions and not artifact.is_latest:
        delete_in_storage = False
        if storage:
            logger.warning(
                "storage argument is ignored; can't delete store of a previous version if overwrite_versions is True"
            )
    elif artifact.key is None or (
        artifact._key_is_virtual and artifact._real_key is None
    ):
        # do not ask for confirmation also if storage is None
        delete_in_storage = storage is None or storage
    else:
        # for artifacts with non-virtual semantic storage keys (key is not None)
        # ask for extra-confirmation if storage is None
        # the wording here is critical to avoid accidental deletions
        if storage is None:
            response = input(
                f"Artifact record deleted. Do you ALSO want to delete the data in storage at {path}? (y/n) You can't undo"
                " this action."
            )
            delete_in_storage = response == "y"
        else:
            delete_in_storage = storage
    if not delete_in_storage:
        logger.important(f"a file/folder remains here: {path}")
    # we don't yet have logic to bring back the deleted metadata record
    # in case storage deletion fails - this is important for ACID down the road
    if delete_in_storage:
        delete_msg = delete_storage(path, raise_file_not_found_error=False)
        if delete_msg != "did-not-delete":
            logger.success(f"deleted {colors.yellow(f'{path}')}")


class LazyArtifact:
    """Lazy artifact for streaming to auto-generated internal paths.

    This is needed when it is desirable to stream to a `lamindb` auto-generated internal path
    and register the path as an artifact (see :class:`~lamindb.Artifact`).

    This object creates a real artifact on `.save()` with the provided arguments.

    Args:
        suffix: The suffix for the auto-generated internal path
        overwrite_versions: Whether to overwrite versions.
        **kwargs: Keyword arguments for the artifact to be created.

    Examples:

        Create a lazy artifact, write to the path and save to get a real artifact::

            lazy = ln.Artifact.from_lazy(suffix=".zarr", overwrite_versions=True, key="mydata.zarr")
            zarr.open(lazy.path, mode="w")["test"] = np.array(["test"]) # stream to the path
            artifact = lazy.save()
    """

    def __init__(self, suffix: str, overwrite_versions: bool, **kwargs):
        self.kwargs = kwargs
        self.kwargs["overwrite_versions"] = overwrite_versions

        if (key := kwargs.get("key")) is not None and extract_suffix_from_path(
            PurePosixPath(key)
        ) != suffix:
            raise ValueError(
                "The suffix argument and the suffix of key should be the same."
            )

        uid, _ = create_uid(n_full_id=20)
        storage_key = auto_storage_key_from_artifact_uid(
            uid, suffix, overwrite_versions=overwrite_versions
        )
        storepath = setup_settings.storage.root / storage_key

        self._path = storepath

    @property
    def path(self) -> UPath:
        return self._path

    def save(self, upload: bool | None = None, **kwargs) -> Artifact:
        artifact = Artifact(self.path, _is_internal_call=True, **self.kwargs)
        return artifact.save(upload=upload, **kwargs)

    def __repr__(self) -> str:  # pragma: no cover
        show_kwargs = {k: v for k, v in self.kwargs.items() if v is not None}
        return (
            f"LazyArtifact object with\n path: {self.path}\n arguments: {show_kwargs}"
        )


class Artifact(SQLRecord, IsVersioned, TracksRun, TracksUpdates):
    """Datasets & models stored as files, folders, or arrays.

    Some artifacts are table- or array-like, e.g., when stored as `.parquet`, `.h5ad`, `.zarr`, or `.tiledb`.

    Args:
        path: `UPathStr` A path to a local or remote folder or file from which to create the artifact.
        key: `str | None = None` A key within the storage location, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a version family.
        description: `str | None = None` A description.
        kind: `Literal["dataset", "model"] | str | None = None` Distinguish models from datasets from other files & folders.
        features: `dict | None = None` External features to annotate the artifact with via :class:`~lamindb.models.FeatureManager.set_values`.
        schema: `Schema | None = None` A schema to validate features.
        revises: `Artifact | None = None` Previous version of the artifact. An alternative to passing `key` when creating a new version.
        overwrite_versions: `bool | None = None` Whether to overwrite versions. Defaults to `True` for folders and `False` for files.
        run: `Run | bool | None = None` The run that creates the artifact. If `False`, suppress tracking the run.
            If `None`, infer the run from the global run context.
        branch: `Branch | None = None` The branch of the artifact. If `None`, uses the current branch.
        space: `Space | None = None` The space of the artifact. If `None`, uses the current space.
        storage: `Storage | None = None` The storage location for the artifact. If `None`, uses the default storage location.
            You can see and set the default storage location in :attr:`~lamindb.core.Settings.storage`.
        skip_hash_lookup: Skip the hash lookup so that a new artifact is created even if an identical artifact already exists.

    Examples:

        Create an artifact **from a local file or folder**::

            artifact = ln.Artifact("./my_file.parquet", key="examples/my_file.parquet").save()
            artifact = ln.Artifact("./my_folder", key="project1/my_folder").save()

        Calling `.save()` copies or uploads the file to the default storage location of your lamindb instance.
        If you create an artifact **from a remote file or folder**, lamindb registers the S3 `key` and avoids copying the data::

            artifact = ln.Artifact("s3://my_bucket/my_folder/my_file.csv").save()

        If you then want to query & access the artifact later on, this is how you do it::

            artifact = ln.Artifact.get(key="examples/my_file.parquet")
            cached_path = artifact.cache()  # sync to local cache & get local path

        If the storage format supports it, you can load the artifact directly into memory or query it through a streaming interface, e.g., for parquet files::

            df = artifact.load()               # load parquet file as DataFrame
            pyarrow_dataset = artifact.open()  # open a streaming file-like object

        If you want to **validate & annotate** a dataframe or an array using the feature & label registries,
        pass `schema` to one of the `.from_dataframe()`, `.from_anndata()`, ... constructors::

            artifact = ln.Artifact.from_dataframe(
                "./my_file.parquet",
                key="my_dataset.parquet",
                schema="valid_features"
            ).save()

        To annotate by **external features**::

            artifact = ln.Artifact("./my_file.parquet", features={"cell_type_by_model": "T cell"}).save()

        You can make a **new version** of an artifact by passing an existing `key`::

            artifact_v2 = ln.Artifact("./my_file.parquet", key="examples/my_file.parquet").save()
            artifact_v2.versions.to_dataframe()  # see all versions

        You can write artifacts to **non-default storage locations** by passing the `storage` argument::

            storage_loc = ln.Storage.get(root="s3://my_bucket")  # get storage location, or create via ln.Storage(root="s3://my_bucket").save()
            ln.Artifact("./my_file.parquet", key="examples/my_file.parquet", storage=storage_loc).save()  # upload to s3://my_bucket

        Sometimes you want to **avoid mapping the artifact into a path hierarchy**, and you only pass `description`::

            artifact = ln.Artifact("./my_folder", description="My folder").save()
            artifact_v2 = ln.Artifact("./my_folder", revises=old_artifact).save()  # need to version based on `revises`, a shared description does not trigger a new version

    Notes:

        .. dropdown:: Typical storage formats & their API accessors

            Arrays:

            - Table: `.csv`, `.tsv`, `.parquet`, `.ipc` ⟷ `DataFrame`, `pyarrow.Table`
            - Annotated matrix: `.h5ad`, `.h5mu`, `.zrad` ⟷ `AnnData`, `MuData`
            - Generic array: HDF5 group, zarr group, TileDB store ⟷ HDF5, zarr, TileDB loaders

            Non-arrays:

            - Image: `.jpg`, `.png` ⟷ `np.ndarray`, ...
            - Fastq: `.fastq` ⟷ /
            - VCF: `.vcf` ⟷ /
            - QC: `.html` ⟷ /

            You'll find these values in the `suffix` & `otype` (object type) fields.

            LaminDB makes some default choices (e.g., serialize a `DataFrame` as a `.parquet` file).

        .. dropdown:: Will artifacts get duplicated?

            If an artifact with the exact same hash already exists, `Artifact()` returns the existing artifact.

            In concurrent workloads where the same artifact is created repeatedly at the exact same time, `.save()`
            detects the duplication and will return the existing artifact.

        .. dropdown:: Why does the constructor look the way it looks?

            It's inspired by APIs building on AWS S3.

            Both boto3 and quilt select a bucket (a storage location in LaminDB) and define a target path through a `key` argument.

            In `boto3 <https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/s3/bucket/upload_file.html>`__::

                # signature: S3.Bucket.upload_file(filepath, key)
                import boto3
                s3 = boto3.resource('s3')
                bucket = s3.Bucket('mybucket')
                bucket.upload_file('/tmp/hello.txt', 'hello.txt')

            In `quilt3 <https://docs.quiltdata.com/api-reference/bucket>`__::

                # signature: quilt3.Bucket.put_file(key, filepath)
                import quilt3
                bucket = quilt3.Bucket('mybucket')
                bucket.put_file('hello.txt', '/tmp/hello.txt')

    See Also:
        :class:`~lamindb.Storage`
            Storage locations for artifacts.
        :class:`~lamindb.Collection`
            Collections of artifacts.
        :meth:`~lamindb.Artifact.from_dataframe`
            Create an artifact from a `DataFrame`.
        :meth:`~lamindb.Artifact.from_anndata`
            Create an artifact from an `AnnData`.

    """

    class Meta(SQLRecord.Meta, IsVersioned.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False
        app_label = "lamindb"
        constraints = [
            # a simple hard unique constraint on `hash` clashes with the fact
            # that pipelines sometimes aim to ingest the exact same file in different
            # folders
            # the conditional composite constraint allows duplicating files in different parts of the
            # file hierarchy, but errors if the same file is to be registered with the same key
            # In SQL, NULL values are treated specially in unique constraints.
            # Multiple NULL values are not considered equal to each other for uniqueness purposes.
            # For non-NULL keys
            models.UniqueConstraint(
                fields=["storage", "key", "hash"],
                condition=models.Q(key__isnull=False),
                name="unique_artifact_storage_key_hash_not_null",
            ),
            # For NULL keys (only storage + hash need to be unique)
            models.UniqueConstraint(
                fields=["storage", "hash"],
                condition=models.Q(key__isnull=True),
                name="unique_artifact_storage_hash_null_key",
            ),
        ]

    _aux_fields: dict[str, tuple[str, type]] = {
        "0": ("_is_saved_to_storage_location", bool),
    }
    _len_full_uid: int = 20
    _len_stem_uid: int = 16
    _name_field: str = "key"

    @property
    def features(self) -> FeatureManager:
        """Feature manager.

        Typically, you annotate a dataset with features by defining a `Schema` and passing it to the `Artifact` constructor.

        Here is how to do annotate an artifact ad hoc::

            artifact.features.add_values({
                "species": organism,  # here, organism is an Organism record
                "scientist": ['Barbara McClintock', 'Edgar Anderson'],
                "temperature": 27.6,
                "experiment": "Experiment 1"
            })

        Query artifacts by features::

            ln.Artifact.filter(scientist="Barbara McClintock")

        Features may or may not be part of the dataset, i.e., the artifact content in storage.
        For instance, the :class:`~lamindb.curators.DataFrameCurator` flow validates the columns of a
        `DataFrame`-like artifact and annotates it with features corresponding to these columns.
        `artifact.features.add_values`, by contrast, does not validate the content of the artifact.

        .. dropdown:: An example for a model-like artifact

            ::

                artifact.features.add_values({
                    "hidden_size": 32,
                    "bottleneck_size": 16,
                    "batch_size": 32,
                    "preprocess_params": {
                        "normalization_type": "cool",
                        "subset_highlyvariable": True,
                    },
                })

        To validate external features::

            schema = ln.Schema([ln.Feature(name="species", dtype=str).save()]).save()
            artifact.features.add_values({"species": "bird"}, schema=schema)
        """
        from ._feature_manager import FeatureManager

        return FeatureManager(self)

    @property
    def labels(self) -> LabelManager:
        """Label manager.

        A way to access all label annotations of an artifact, irrespective of their type.

        To annotate with labels, use the type-specific accessor,
        for example::

            experiment = ln.Record(name="Experiment 1").save()
            artifact.records.add(experiment)
            project = ln.Project(name="Project A").save()
            artifact.projects.add(project)
        """
        from ._label_manager import LabelManager

        return LabelManager(self)

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=_len_full_uid
    )
    """A universal random id."""
    # the max length of 1024 equals the max length of a S3 key
    key: str | None = CharField(db_index=True, null=True, max_length=1024)
    """A (virtual) relative file path within the artifact's storage location.

    Setting a `key` is useful to automatically group artifacts into a version family.

    LaminDB defaults to a virtual file path to make renaming of data in object storage easy.

    If you register existing files in a storage location, the `key` equals the
    actual filepath on the underyling filesytem or object store.
    """
    _real_key: str | None = CharField(db_index=True, null=True, max_length=1024)
    """An optional real storage key."""
    # db_index on description because sometimes we query for equality in the case of artifacts
    description: str | None = TextField(null=True, db_index=True)
    """A description."""
    storage: Storage = ForeignKey(
        Storage, PROTECT, related_name="artifacts", editable=False
    )
    """Storage location, e.g. an S3 or GCP bucket or a local directory."""
    suffix: str = CharField(max_length=30, db_index=True, editable=False)
    # Initially, we thought about having this be nullable to indicate folders
    # But, for instance, .zarr is stored in a folder that ends with a .zarr suffix
    """Path suffix or empty string if no canonical suffix exists.

    This is either a file suffix (`".csv"`, `".h5ad"`, etc.) or the empty string "".
    """
    kind: ArtifactKind | str | None = CharField(
        max_length=20,
        db_index=True,
        null=True,
    )
    """:class:`~lamindb.base.types.ArtifactKind` or custom `str` value (default `None`)."""
    otype: str | None = CharField(
        max_length=64, db_index=True, null=True, editable=False
    )
    """Default Python object type, e.g., DataFrame, AnnData."""
    size: int | None = BigIntegerField(
        null=True, db_index=True, default=None, editable=False
    )
    """Size in bytes.

    Examples: 1KB is 1e3 bytes, 1MB is 1e6, 1GB is 1e9, 1TB is 1e12 etc.
    """
    hash: str | None = CharField(
        max_length=HASH_LENGTH, db_index=True, null=True, editable=False
    )
    """Hash or pseudo-hash of artifact content.

    Useful to ascertain integrity and avoid duplication.
    """
    n_files: int | None = BigIntegerField(
        null=True, db_index=True, default=None, editable=False
    )
    """Number of files for folder-like artifacts, `None` for file-like artifacts.

    Note that some arrays are also stored as folders, e.g., `.zarr` or `.tiledbsoma`.
    """
    n_observations: int | None = BigIntegerField(
        null=True, db_index=True, default=None, editable=False
    )
    """Number of observations.

    Typically, this denotes the first array dimension.
    """
    _hash_type: str | None = CharField(
        max_length=30, db_index=True, null=True, editable=False
    )
    """Type of hash."""
    run: Run | None = ForeignKey(
        Run,
        PROTECT,
        related_name="output_artifacts",
        null=True,
        default=None,
        editable=False,
    )
    """Run that created the artifact."""
    input_of_runs: Run = models.ManyToManyField(Run, related_name="input_artifacts")
    """Runs that use this artifact as an input."""
    _subsequent_runs: Run = models.ManyToManyField(
        "Run",
        related_name="_recreated_artifacts",
        db_table="lamindb_artifact__previous_runs",  # legacy name, change in lamindb v2
    )
    """Runs that re-created the record after initial creation."""
    collections: Collection
    """The collections that this artifact is part of."""
    schema: Schema | None = ForeignKey(
        Schema,
        PROTECT,
        null=True,
        default=None,
        related_name="validated_artifacts",
    )
    """The schema that validated this artifact in a :class:`~lamindb.curators.core.Curator`."""
    feature_sets: Schema = models.ManyToManyField(
        Schema, related_name="artifacts", through="ArtifactSchema"
    )
    """The feature sets measured by the artifact."""
    _feature_values: FeatureValue = models.ManyToManyField(
        FeatureValue, through="ArtifactFeatureValue", related_name="artifacts"
    )
    """Non-categorical feature values for annotation."""
    _key_is_virtual: bool = BooleanField()
    """Indicates whether `key` is virtual or part of an actual file path."""
    # be mindful that below, passing related_name="+" leads to errors
    _actions: Artifact = models.ManyToManyField(
        "self", symmetrical=False, related_name="_action_targets"
    )
    """Actions to attach for the UI."""
    created_by: User = ForeignKey(
        "lamindb.User",
        PROTECT,
        default=current_user_id,
        related_name="created_artifacts",
        editable=False,
    )
    """Creator of record."""
    _overwrite_versions: bool = BooleanField(default=None)
    """See corresponding property `overwrite_versions`."""
    ulabels: ULabel
    """The ulabels annotating this artifact."""
    users: User
    """The users annotating this artifact."""
    projects: Project
    """The projects annotating this artifact."""
    references: Reference
    """The references annotating this artifact."""
    records: Record
    """The records annotating this artifact."""
    linked_in_records: Record = models.ManyToManyField(
        "Record", through="RecordArtifact", related_name="linked_artifacts"
    )
    """This artifact is linked in these records as a value."""
    blocks: ArtifactBlock
    """The blocks that annotate this artifact."""

    @overload
    def __init__(
        self,
        # we're not choosing the name "path" for this arg because
        # it could be confused with `artifact.path`
        # "data" conveys better that this is input data that's ingested
        # and will be moved to a target path at `artifact.path`
        data: UPathStr,
        *,
        key: str | None = None,
        description: str | None = None,
        kind: ArtifactKind | str | None = None,
        features: dict[str, Any] | None = None,
        schema: Schema | None = None,
        revises: Artifact | None = None,
        overwrite_versions: bool | None = None,
        run: Run | False | None = None,
        storage: Storage | None = None,
        branch: Branch | None = None,
        space: Space | None = None,
        skip_hash_lookup: bool = False,
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
        # check whether we are called with db args
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        # now proceed with the user-facing constructor
        if len(args) > 1:
            raise ValueError("Only one non-keyword arg allowed: path")

        if "data" in kwargs:
            warnings.warn(
                "`data` argument was renamed to `path` and will be removed in a future release.",
                DeprecationWarning,
                stacklevel=2,
            )
            path = kwargs.pop("data")
        else:
            path = kwargs.pop("path") if len(args) == 0 else args[0]
        kind: str = kwargs.pop("kind", None)
        key: str | None = kwargs.pop("key", None)
        run_id: int | None = kwargs.pop("run_id", None)  # for REST API
        run: Run | None = kwargs.pop("run", None)
        using_key = kwargs.pop("using_key", None)
        description: str | None = kwargs.pop("description", None)
        revises: Artifact | None = kwargs.pop("revises", None)
        overwrite_versions: bool | None = kwargs.pop("overwrite_versions", None)
        version: str | None = kwargs.pop("version", None)
        schema: Schema | None = kwargs.pop("schema", None)
        features: dict[str, Any] | None = kwargs.pop("features", None)
        skip_hash_lookup: bool = kwargs.pop("skip_hash_lookup", False)
        to_disk_kwargs: dict[str, Any] | None = kwargs.pop("to_disk_kwargs", None)

        # validate external features if passed with a schema
        if features is not None:
            self._external_features = features
            if schema is not None:
                from lamindb.curators.core import ExperimentalDictCurator

                validation_schema = schema
                ExperimentalDictCurator(features, validation_schema).validate()

        branch = kwargs.pop("branch", None)
        assert "branch_id" not in kwargs, "Please pass branch instead of branch_id."  # noqa: S101
        space = kwargs.pop("space", None)
        assert "space_id" not in kwargs, "Please pass space instead of space_id."  # noqa: S101
        format = kwargs.pop("format", None)
        _is_internal_call = kwargs.pop("_is_internal_call", False)
        skip_check_exists = kwargs.pop("skip_check_exists", False)
        storage_was_passed = False
        if "storage" in kwargs:
            storage = kwargs.pop("storage")
            storage_was_passed = True
        elif (
            setup_settings.instance.keep_artifacts_local
            and setup_settings.instance._local_storage is not None
        ):
            storage = setup_settings.instance.local_storage.record
        else:
            storage = setup_settings.instance.storage.record
        if space is None:
            from lamindb import context as run_context

            if run_context.space is not None:
                space = run_context.space
            elif setup_settings.space is not None:
                space = setup_settings.space
        if space is not None and space != storage.space:
            if storage_was_passed:
                logger.warning(
                    "storage argument ignored as storage information from space takes precedence"
                )
            storage_locs_for_space = Storage.filter(space=space)
            n_storage_locs_for_space = len(storage_locs_for_space)
            if n_storage_locs_for_space == 0:
                raise NoStorageLocationForSpace(
                    "No storage location found for space.\n"
                    "Either create one via ln.Storage(root='create-s3', space=space).save()\n"
                    "Or start managing access to an existing storage location via the space: storage_loc.space = space; storage.save()"
                )
            else:
                storage = storage_locs_for_space.first()
                if n_storage_locs_for_space > 1:
                    logger.warning(
                        f"more than one storage location for space {space}, choosing {storage}"
                    )
        otype = kwargs.pop("otype") if "otype" in kwargs else None
        if isinstance(path, str) and path.startswith("s3:///"):
            # issue in Groovy / nf-lamin producing malformed S3 paths
            # https://laminlabs.slack.com/archives/C08J590666Q/p1751315027830849?thread_ts=1751039961.479259&cid=C08J590666Q
            path = path.replace("s3:///", "s3://")
        otype = _check_otype_artifact(
            data=path, otype=otype, cloud_warning=not _is_internal_call
        )
        if "type" in kwargs:
            logger.warning("`type` will be removed soon, please use `kind`")
            kind = kwargs.pop("type")
        if not len(kwargs) == 0:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Artifact)])
            raise FieldValidationError(
                f"Only {valid_keywords} can be passed, you passed: {kwargs}"
            )
        if revises is not None and key is not None and revises.key != key:
            logger.warning(f"renaming artifact from '{revises.key}' to {key}")
        if revises is not None:
            if not isinstance(revises, Artifact):
                raise TypeError("`revises` has to be of type `Artifact`")
            if description is None:
                description = revises.description
        if key is not None and AUTO_KEY_PREFIX in key:
            raise ValueError(
                f"Do not pass key that contains a managed storage path in `{AUTO_KEY_PREFIX}`"
            )
        # below is for internal calls that require defining the storage location
        # ahead of constructing the Artifact
        if isinstance(path, (str, Path)) and AUTO_KEY_PREFIX in str(path):
            if _is_internal_call:
                is_automanaged_path = True
                user_provided_key = key
                key = None
            else:
                raise ValueError(
                    f"Do not pass path inside the `{AUTO_KEY_PREFIX}` directory."
                )
        else:
            is_automanaged_path = False

        provisional_uid, revises = create_uid(revises=revises, version=version)
        run = get_run(run)
        kwargs_or_artifact, privates = get_artifact_kwargs_from_data(
            data=path,
            key=key,
            run=run,
            format=format,
            provisional_uid=provisional_uid,
            version=version,
            storage=storage,
            using_key=using_key,
            skip_check_exists=skip_check_exists,
            overwrite_versions=overwrite_versions,
            skip_hash_lookup=skip_hash_lookup,
            to_disk_kwargs=to_disk_kwargs,
        )

        # an object with the same hash already exists
        if isinstance(kwargs_or_artifact, Artifact):
            from .sqlrecord import init_self_from_db, update_attributes

            init_self_from_db(self, kwargs_or_artifact)
            # update key from inferred value
            key = privates
            # adding "key" here is dangerous because key might be auto-populated
            attr_to_update = {"description": description}
            if schema is not None:
                attr_to_update["schema"] = schema
            if kwargs_or_artifact._key_is_virtual and kwargs_or_artifact.key is None:
                attr_to_update["key"] = key
            elif self.key != key and key is not None:
                if not self.path.exists():
                    logger.warning(f"updating previous key {self.key} to new key {key}")
                    self.key = key
                    assert self.path.exists(), (  # noqa: S101
                        f"The underlying file for artifact {self} does not exist anymore, clean up the artifact record."
                    )  # noqa: S101
                    self._skip_key_change_check = (
                        True  # otherwise not allowed to change real keys
                    )
                else:
                    logger.warning(
                        f"key {self.key} on existing artifact differs from passed key {key}, keeping original key; update manually if needed or pass skip_hash_lookup if you want to duplicate the artifact"
                    )
            update_attributes(self, attr_to_update)
            populate_subsequent_run(self, run)
            return None
        else:
            kwargs = kwargs_or_artifact
            kwargs["schema"] = schema

        if revises is None:
            revises = kwargs_or_artifact.pop("revises")

        if path is not None:
            self._local_filepath = privates["local_filepath"]
            self._cloud_filepath = privates["cloud_filepath"]
            self._memory_rep = privates["memory_rep"]
            self._to_store = not privates["check_path_in_storage"]

        if is_automanaged_path and _is_internal_call:
            kwargs["_key_is_virtual"] = True
            assert AUTO_KEY_PREFIX in kwargs["key"]  # noqa: S101
            uid = (
                kwargs["key"].replace(AUTO_KEY_PREFIX, "").replace(kwargs["suffix"], "")
            )
            kwargs["key"] = user_provided_key
            if revises is not None:
                assert uid.startswith(revises.stem_uid)  # noqa: S101
            if len(uid) == 16:
                if revises is None:
                    uid += "0000"
                else:
                    uid, revises = create_uid(revises=revises, version=version)
            kwargs["uid"] = uid

        # only set key now so that we don't perform a look-up on it in case revises is passed
        if revises is not None and revises.key is not None and kwargs["key"] is None:
            kwargs["key"] = revises.key

        if run_id is not None:
            kwargs["run_id"] = run_id
        kwargs["kind"] = kind
        kwargs["version"] = version
        kwargs["description"] = description
        kwargs["branch"] = branch
        kwargs["space"] = space
        kwargs["otype"] = otype
        kwargs["revises"] = revises
        # this check needs to come down here because key might be populated from an
        # existing file path during get_artifact_kwargs_from_data()
        if (
            kwargs["key"] is None
            and kwargs["description"] is None
            and kwargs["run"] is None
        ):
            raise ValueError("Pass one of key, run or description as a parameter")

        super().__init__(**kwargs)

    @property
    @deprecated("kind")
    def type(self) -> str:
        return self.kind

    @property
    @deprecated("otype")
    def _accessor(self) -> str:
        return self.otype

    @property
    @deprecated("features")
    def params(self) -> str:
        return self.features

    @property
    def transform(self) -> Transform | None:
        """Transform whose run created the artifact."""
        return self.run.transform if self.run is not None else None

    @property
    @deprecated("n_files")
    def n_objects(self) -> int:
        return self.n_files

    @property
    def overwrite_versions(self) -> bool:
        """Indicates whether to keep or overwrite versions.

        It defaults to `False` for file-like artifacts and to `True` for folder-like artifacts.

        Note that this requires significant storage space for large folders with
        many duplicated files. Currently, `lamindb` does *not* de-duplicate files across
        versions as in git, but keeps all files for all versions of the folder in storage.
        """
        return self._overwrite_versions

    @property
    def path(self) -> Path:
        """Path.

        Example::

            import lamindb as ln

            # File in cloud storage, here AWS S3:
            artifact = ln.Artifact("s3://my-bucket/my-file.csv").save()
            artifact.path
            #S3QueryPath('s3://my-bucket/my-file.csv')

            # File in local storage:
            ln.Artifact("./myfile.csv", key="myfile.csv").save()
            artifact.path
            #> PosixPath('/home/runner/work/lamindb/lamindb/docs/guide/mydata/myfile.csv')
        """
        filepath, _ = filepath_from_artifact(self, using_key=settings._using_key)
        return filepath

    @property
    def _cache_path(self) -> UPath:
        filepath, cache_key = filepath_cache_key_from_artifact(
            self, using_key=settings._using_key
        )
        if isinstance(filepath, LocalPathClasses):
            return filepath
        return setup_settings.paths.cloud_to_local_no_update(
            filepath, cache_key=cache_key
        )

    @classmethod
    def get(
        cls,
        idlike: int | str | None = None,
        *,
        key: str | None = None,
        path: str | None = None,
        is_run_input: bool | Run = False,
        **expressions,
    ) -> Artifact:
        """Get a single artifact.

        Args:
            idlike: Either a uid stub, uid or an integer id.
            key: An optional key to query for.
            path: An optional full path to query for, including the storage root.
            is_run_input: Whether to track this artifact as run input.
            expressions: Other fields and values passed as Django query expressions.

        Raises:
            :exc:`lamindb.errors.DoesNotExist`: In case no matching record is found.

        See Also:
            - Guide: :doc:`registries`
            - Method in `SQLRecord` base class: :meth:`~lamindb.models.SQLRecord.get`

        Examples:

            ::

                artifact = ln.Artifact.get("tCUkRcaEjTjhtozp")       # gets latest version for family tCUkRcaEjTjhtozp
                artifact = ln.Artifact.get("tCUkRcaEjTjhtozp0005")   # gets version 0005 for family tCUkRcaEjTjhtozp
                artifact = ln.Artifact.get(key="examples/my_file.parquet")               # gets latest version for a key
                artifact = ln.Artifact.get(key="examples/my_file.parquet", version="2")  # pass a version tag
                artifact = ln.Artifact.get(path="s3://bucket/folder/adata.h5ad")
        """
        if key is not None:
            expressions["key"] = key
        if path is not None:
            expressions["path"] = path
        return QuerySet(model=cls).get(idlike, is_run_input=is_run_input, **expressions)

    @classmethod
    def filter(
        cls,
        *queries,
        **expressions,
    ) -> QuerySet:
        """Query a set of artifacts.

        Args:
            *queries: `Q` expressions.
            **expressions: Features & fields via the Django query syntax.

        See Also:
            - Guide: :doc:`docs:registries`

        Examples:

            Query by fields::

                ln.Arfifact.filter(key="examples/my_file.parquet")

            Query by features::

                ln.Arfifact.filter(cell_type_by_model__name="T cell")

        """
        # from Registry metaclass
        return type(cls).filter(cls, *queries, **expressions)

    @classmethod
    def from_lazy(
        cls,
        suffix: str,
        overwrite_versions: bool,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        **kwargs,
    ) -> LazyArtifact:
        """Create a lazy artifact for streaming to auto-generated internal paths.

        This is needed when it is desirable to stream to a `lamindb` auto-generated internal path
        and register the path as an artifact.

        The lazy artifact object (see :class:`~lamindb.models.LazyArtifact`) creates a real artifact
        on `.save()` with the provided arguments.

        Args:
            suffix: The suffix for the auto-generated internal path
            overwrite_versions: Whether to overwrite versions.
            key: An optional key to reference the artifact.
            description: A description.
            run: The run that creates the artifact.
            **kwargs: Other keyword arguments for the artifact to be created.

        Examples:

            Create a lazy artifact, write to the path and save to get a real artifact::

                lazy = ln.Artifact.from_lazy(suffix=".zarr", overwrite_versions=True, key="mydata.zarr")
                zarr.open(lazy.path, mode="w")["test"] = np.array(["test"]) # stream to the path
                artifact = lazy.save()
        """
        args = {"key": key, "description": description, "run": run, **kwargs}
        return LazyArtifact(suffix, overwrite_versions, **args)

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        *,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        schema: Schema | Literal["valid_features"] | None = None,
        features: dict[str, Any] | None = None,
        parquet_kwargs: dict[str, Any] | None = None,
        csv_kwargs: dict[str, Any] | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from `DataFrame`, optionally validate & annotate.

        Args:
            df: A `DataFrame` object.
            key: A relative path within default storage,
                e.g., `"myfolder/myfile.parquet"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.
            schema: A schema that defines how to validate & annotate.
            features: Additional external features to annotate the artifact via :class:`~lamindb.models.FeatureManager.set_values`.
            parquet_kwargs: Additional keyword arguments passed to the
                `pandas.DataFrame.to_parquet` method, which are passed
                on to `pyarrow.parquet.ParquetWriter`.
            csv_kwargs: Additional keyword arguments passed to the
                `pandas.DataFrame.to_csv` method.

        Examples:

            No validation and annotation::

                ln.Artifact.from_dataframe(df, key="examples/dataset1.parquet").save()

            With validation and annotation::

                ln.Artifact.from_dataframe(df, key="examples/dataset1.parquet", schema="valid_features").save()

            Under-the-hood, this uses the following build-in schema (:func:`~lamindb.examples.schemas.valid_features`)::

                schema = ln.Schema(name="valid_features", itype="Feature").save()

            External features:

            .. literalinclude:: scripts/curate_dataframe_external_features.py
               :language: python

            Parquet kwargs:

            .. literalinclude:: scripts/test_artifact_parquet.py
               :language: python
        """
        from lamindb import examples

        if "format" not in kwargs and key is not None and key.endswith(".csv"):
            kwargs["format"] = ".csv"
        if schema == "valid_features":
            schema = examples.schemas.valid_features()
        to_disk_kwargs: dict[str, Any] = parquet_kwargs or csv_kwargs
        artifact = Artifact(  # type: ignore
            data=df,
            key=key,
            run=run,
            description=description,
            revises=revises,
            otype="DataFrame",
            kind="dataset",
            to_disk_kwargs=to_disk_kwargs,
            **kwargs,
        )
        artifact.n_observations = len(df)
        if features is not None:
            artifact._external_features = features
        if schema is not None:
            from lamindb.curators.core import DataFrameCurator

            if not artifact._state.adding and artifact.suffix != ".parquet":
                logger.warning(
                    f"not re-validating existing artifact as it was stored as {artifact.suffix}, "
                    "which does not maintain categorical dtype information"
                )
                return artifact

            curator = DataFrameCurator(artifact, schema, features=features)
            curator.validate()
            artifact.schema = schema
            artifact._curator = curator
        return artifact

    @classmethod
    @deprecated("from_dataframe")
    def from_df(
        cls,
        df: pd.DataFrame,
        *,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        schema: Schema | None = None,
        **kwargs,
    ) -> Artifact:
        return cls.from_dataframe(
            df,
            key=key,
            description=description,
            run=run,
            revises=revises,
            schema=schema,
            **kwargs,
        )

    @classmethod
    def from_anndata(
        cls,
        adata: Union[AnnData, UPathStr],
        *,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        schema: Schema
        | Literal["ensembl_gene_ids_and_valid_features_in_obs"]
        | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from `AnnData`, optionally validate & annotate.

        Args:
            adata: An `AnnData` object or a path of AnnData-like.
            key: A relative path within default storage,
                e.g., `"myfolder/myfile.h5ad"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.
            schema: A schema that defines how to validate & annotate.

        See Also:
            :meth:`~lamindb.Collection`
                Track collections.
            :class:`~lamindb.Feature`
                Track features.

        Example:

            No validation and annotation::

                ln.Artifact.from_anndata(adata, key="examples/dataset1.h5ad").save()

            With validation and annotation::

                ln.Artifact.from_dataframe(adata, key="examples/dataset1.h5ad", schema="ensembl_gene_ids_and_valid_features_in_obs").save()

            Under-the-hood, this uses the following build-in schema (:func:`~lamindb.examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs`):

            .. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py
               :language: python

            This schema tranposes the `var` DataFrame during curation, so that one validates and annotates the columns of `var.T`, i.e., `[ENSG00000153563, ENSG00000010610, ENSG00000170458]`.
            If one doesn't transpose, one would annotate the columns of `var`, i.e., `[gene_symbol, gene_type]`.

            .. image:: https://lamin-site-assets.s3.amazonaws.com/.lamindb/gLyfToATM7WUzkWW0001.png
               :width: 800px

        """
        from lamindb import examples

        if not data_is_scversedatastructure(adata, "AnnData"):
            raise ValueError(
                "data has to be an AnnData object or a path to AnnData-like"
            )
        if schema == "ensembl_gene_ids_and_valid_features_in_obs":
            schema = (
                examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs()
            )
        _anndata_n_observations(adata)
        artifact = Artifact(  # type: ignore
            data=adata,
            key=key,
            run=run,
            description=description,
            revises=revises,
            otype="AnnData",
            kind="dataset",
            **kwargs,
        )
        # this is done instead of _anndata_n_observations(adata)
        # because we need a proper path through create_path for cloud paths
        # for additional upath options etc that create_path adds
        obj_for_obs: AnnData | UPath
        if hasattr(artifact, "_memory_rep") and artifact._memory_rep is not None:
            obj_for_obs = artifact._memory_rep
        else:
            # returns ._local_filepath for local files
            # and the proper path through create_path for cloud paths
            obj_for_obs = artifact.path
        artifact.n_observations = _anndata_n_observations(obj_for_obs)
        if schema is not None:
            from ..curators import AnnDataCurator

            curator = AnnDataCurator(artifact, schema)
            curator.validate()
            artifact.schema = schema
            artifact._curator = curator
        return artifact

    @classmethod
    def from_mudata(
        cls,
        mdata: Union[MuData, UPathStr],
        *,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        schema: Schema | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from `MuData`, optionally validate & annotate.

        Args:
            mdata: A `MuData` object.
            key: A relative path within default storage,
                e.g., `"myfolder/myfile.h5mu"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.
            schema: A schema that defines how to validate & annotate.

        See Also:
            :meth:`~lamindb.Collection`
                Track collections.
            :class:`~lamindb.Feature`
                Track features.

        Example::

            import lamindb as ln

            mdata = ln.examples.datasets.mudata_papalexi21_subset()
            artifact = ln.Artifact.from_mudata(mdata, key="mudata_papalexi21_subset.h5mu").save()
        """
        if not data_is_scversedatastructure(mdata, "MuData"):
            raise ValueError("data has to be a MuData object or a path to MuData-like")
        artifact = Artifact(  # type: ignore
            data=mdata,
            key=key,
            run=run,
            description=description,
            revises=revises,
            otype="MuData",
            kind="dataset",
            **kwargs,
        )
        if not isinstance(mdata, UPathStr):
            artifact.n_observations = mdata.n_obs
        if schema is not None:
            from ..curators import MuDataCurator

            curator = MuDataCurator(artifact, schema)
            curator.validate()
            artifact.schema = schema
            artifact._curator = curator
        return artifact

    @classmethod
    def from_spatialdata(
        cls,
        sdata: SpatialData | UPathStr,
        *,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        schema: Schema | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from `SpatialData`, optionally validate & annotate.

        Args:
            sdata: A `SpatialData` object.
            key: A relative path within default storage,
                e.g., `"myfolder/myfile.zarr"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.
            schema: A schema that defines how to validate & annotate.

        See Also:
            :meth:`~lamindb.Collection`
                Track collections.
            :class:`~lamindb.Feature`
                Track features.

        Example:

            No validation and annotation::

                import lamindb as ln

                artifact = ln.Artifact.from_spatialdata(sdata, key="my_dataset.zarr").save()

            With validation and annotation.

            .. literalinclude:: scripts/define_schema_spatialdata.py
                :language: python

            .. literalinclude:: scripts/curate_spatialdata.py
                :language: python
        """
        if not data_is_scversedatastructure(sdata, "SpatialData"):
            raise ValueError(
                "data has to be a SpatialData object or a path to SpatialData-like"
            )
        artifact = Artifact(  # type: ignore
            data=sdata,
            key=key,
            run=run,
            description=description,
            revises=revises,
            otype="SpatialData",
            kind="dataset",
            **kwargs,
        )
        # ill-defined https://scverse.zulipchat.com/#narrow/channel/315824-spatial/topic/How.20to.20calculate.20the.20number.20of.20observations.3F
        # artifact.n_observations = ...
        if schema is not None:
            from ..curators import SpatialDataCurator

            curator = SpatialDataCurator(artifact, schema)
            curator.validate()
            artifact.schema = schema
            artifact._curator = curator
        return artifact

    @classmethod
    def from_tiledbsoma(
        cls,
        exp: SOMAExperiment | UPathStr,
        *,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from a `tiledbsoma.Experiment` store.

        Args:
            exp: TileDB-SOMA Experiment object or path to Experiment store.
            key: A relative path within default storage, e.g., `"myfolder/mystore.tiledbsoma"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.

        Example::

            import lamindb as ln

            artifact = ln.Artifact.from_tiledbsoma("s3://mybucket/store.tiledbsoma", description="a tiledbsoma store").save()
        """
        if not data_is_soma_experiment(exp):
            raise ValueError(
                "data has to be a SOMA Experiment object or a path to SOMA Experiment store."
            )

        # SOMAExperiment.uri may have file:// prefix for local paths which needs stripping for filesystem access.
        # Other URI schemes (s3://, etc.) are preserved and supported.
        exp = exp.uri.removeprefix("file://") if not isinstance(exp, UPathStr) else exp

        artifact = Artifact(  # type: ignore
            data=exp,
            key=key,
            run=run,
            description=description,
            revises=revises,
            otype="tiledbsoma",
            kind="dataset",
            **kwargs,
        )
        artifact.n_observations = _soma_n_observations(artifact.path)
        return artifact

    @classmethod
    def from_dir(
        cls,
        path: UPathStr,
        *,
        key: str | None = None,
        run: Run | None = None,
    ) -> SQLRecordList:
        """Create a list of :class:`~lamindb.Artifact` objects from a directory.

        Hint:
            If you have a high number of files (several 100k) and don't want to
            track them individually, create a single :class:`~lamindb.Artifact` via
            ``Artifact(path)`` for them. See, e.g., :doc:`docs:rxrx`.

        Args:
            path: Source path of folder.
            key: Key for storage destination. If `None` and
                directory is in a registered location, the inferred `key` will
                reflect the relative position. If `None` and directory is outside
                of a registered storage location, the inferred key defaults to `path.name`.
            run: A `Run` object.

        Example::

            import lamindb as ln

            dir_path = ln.examples.datasets.generate_cell_ranger_files("sample_001", ln.settings.storage)
            ln.Artifact.from_dir(dir_path).save()  # creates one artifact per file in dir_path
        """
        folderpath: UPath = create_path(path)  # returns Path for local
        storage = settings.storage.record
        using_key = settings._using_key
        storage, use_existing_storage = process_pathlike(folderpath, storage, using_key)
        folder_key_path: PurePath | Path
        if key is None:
            if not use_existing_storage:
                logger.warning(
                    "folder is outside existing storage location, will copy files from"
                    f" {path} to {storage.root}/{folderpath.name}"
                )
                folder_key_path = Path(folderpath.name)
            else:
                # maintain the hierachy within an existing storage location
                folder_key_path = get_relative_path_to_directory(
                    folderpath, UPath(storage.root)
                )
        else:
            folder_key_path = Path(key)

        folder_key = folder_key_path.as_posix()
        # silence fine-grained logging
        verbosity = settings.verbosity
        verbosity_int = settings._verbosity_int
        if verbosity_int >= 1:
            settings.verbosity = "warning"
        artifacts_dict = {}
        for filepath in folderpath.rglob("*"):
            if filepath.is_file():
                relative_path = get_relative_path_to_directory(filepath, folderpath)
                artifact_key = folder_key + "/" + relative_path.as_posix()
                # if creating from rglob, we don't need to check for existence
                artifact = Artifact(
                    filepath, run=run, key=artifact_key, skip_check_exists=True
                )
                artifacts_dict[artifact.uid] = artifact
        settings.verbosity = verbosity

        # run sanity check on hashes
        hashes = [
            artifact.hash
            for artifact in artifacts_dict.values()
            if artifact.hash is not None
        ]
        uids = artifacts_dict.keys()
        n_unique_hashes = len(set(hashes))
        if n_unique_hashes == len(hashes):
            artifacts = SQLRecordList(artifacts_dict.values())
        else:
            # consider exact duplicates (same id, same hash)
            # below can't happen anymore because artifacts is a dict now
            # if len(set(uids)) == len(set(hashes)):
            #     logger.warning("dropping duplicate records in list of artifact records")
            #     artifacts = list(set(uids))
            # consider false duplicates (different id, same hash)
            if not len(set(uids)) == n_unique_hashes:
                seen_hashes = set()
                non_unique_artifacts = {
                    hash: artifact
                    for hash, artifact in artifacts_dict.items()
                    if artifact.hash in seen_hashes or seen_hashes.add(artifact.hash)  # type: ignore
                }
                display_non_unique = "\n    ".join(
                    f"{artifact}" for artifact in non_unique_artifacts
                )
                logger.warning(
                    "there are multiple artifact uids with the same hashes, dropping"
                    f" {len(non_unique_artifacts)} duplicates out of"
                    f" {len(artifacts_dict)} artifacts:\n    {display_non_unique}"
                )
                artifacts = SQLRecordList(
                    [
                        artifact
                        for artifact in artifacts_dict.values()
                        if artifact not in non_unique_artifacts.values()
                    ]
                )
        logger.success(
            f"created {len(artifacts)} artifacts from directory using storage"
            f" {storage.root} and key = {folder_key}/"
        )
        return artifacts

    def replace(
        self,
        data: Union[UPathStr, pd.DataFrame, AnnData, MuData],
        run: Run | None = None,
        format: str | None = None,
    ) -> None:
        """Replace artifact content.

        Args:
            data: A file path.
            run: The run that created the artifact gets
                auto-linked if ``ln.track()`` was called.

        Examples:
            Say we made a change to the content of an artifact, e.g., edited the image
            `paradisi05_laminopathic_nuclei.jpg`.

            This is how we replace the old file in storage with the new file:

            >>> artifact.replace("paradisi05_laminopathic_nuclei.jpg")
            >>> artifact.save()

            Note that this neither changes the storage key nor the filename.

            However, it will update the suffix if it changes.
        """
        storage = settings.storage.record
        run = get_run(run)
        kwargs, privates = get_artifact_kwargs_from_data(
            provisional_uid=self.uid,
            data=data,
            key=self.key,
            run=run,
            format=format,
            storage=storage,
            version=None,
            is_replace=True,
        )

        # this artifact already exists
        if isinstance(kwargs, Artifact):
            return kwargs

        check_path_in_storage = privates["check_path_in_storage"]
        if check_path_in_storage:
            err_msg = (
                "Can only replace with a local path not in any Storage. "
                f"This data is in {Storage.objects.get(id=kwargs['storage_id'])}."
            )
            raise ValueError(err_msg)

        _overwrite_versions = kwargs["_overwrite_versions"]
        if self._overwrite_versions != _overwrite_versions:
            err_msg = "It is not allowed to replace "
            err_msg += "a folder" if self._overwrite_versions else "a file"
            err_msg += " with " + ("a folder." if _overwrite_versions else "a file.")
            raise ValueError(err_msg)

        new_suffix = kwargs["suffix"]
        if new_suffix != self.suffix:
            key = self.key
            real_key = self._real_key
            if key is not None:
                new_key = PurePosixPath(key).with_suffix(new_suffix).as_posix()
            else:
                new_key = None
            if (key is not None and not self._key_is_virtual) or real_key is not None:
                # real_key is not None implies key is not None
                assert key is not None  # noqa: S101
                if real_key is not None:
                    self._clear_storagekey = real_key
                    self._real_key = (
                        PurePosixPath(real_key).with_suffix(new_suffix).as_posix()
                    )
                    warn_msg = f", _real_key '{real_key}' with '{self._real_key}'"
                else:
                    self._clear_storagekey = key
                    warn_msg = ""
                self.key = new_key
                # update old key with the new one so that checks in record pass
                self._old_key = new_key
                logger.warning(
                    f"replacing the file will replace key '{key}' with '{new_key}'{warn_msg}"
                    f" and delete '{self._clear_storagekey}' upon `save()`"
                )
            else:
                # purely virtual key case
                self._clear_storagekey = auto_storage_key_from_artifact(self)
                # might replace None with None, not a big deal
                self.key = new_key
                # update the old key with the new one so that checks in record pass
                self._old_key = new_key

        self.suffix = new_suffix
        self.size = kwargs["size"]
        self.hash = kwargs["hash"]
        self._hash_type = kwargs["_hash_type"]
        self.run_id = kwargs["run_id"]
        self.run = kwargs["run"]
        self.n_files = kwargs["n_files"]

        self._local_filepath = privates["local_filepath"]
        self._cloud_filepath = privates["cloud_filepath"]
        self._memory_rep = privates["memory_rep"]
        # no need to upload if new file is already in storage
        self._to_store = not check_path_in_storage

        # update old suffix with the new one so that checks in record pass
        # replace() supports changing the suffix
        self._old_suffix = self.suffix

    def open(
        self,
        mode: str = "r",
        engine: Literal["pyarrow", "polars"] = "pyarrow",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> (
        AnnDataAccessor
        | SpatialDataAccessor
        | BackedAccessor
        | SOMACollection
        | SOMAExperiment
        | SOMAMeasurement
        | PyArrowDataset
        | Iterator[PolarsLazyFrame]
    ):
        """Open a dataset for streaming.

        Works for `AnnData` (`.h5ad` and `.zarr`), `SpatialData` (`.zarr`),
        generic `hdf5` and `zarr`, `tiledbsoma` objects (`.tiledbsoma`),
        `pyarrow` or `polars` compatible formats
        (`.parquet`, `.csv`, `.ipc` etc. files or directories with such files).

        Args:
            mode: can be `"r"` or `"w"` (write mode) for `tiledbsoma` stores,
                `"r"` or `"r+"` for `AnnData` or `SpatialData` `zarr` stores,
                otherwise should be always `"r"` (read-only mode).
            engine: Which module to use for lazy loading of a dataframe
                from `pyarrow` or `polars` compatible formats.
                This has no effect if the artifact is not a dataframe, i.e.
                if it is an `AnnData,` `hdf5`, `zarr`, `tiledbsoma` object etc.
            is_run_input: Whether to track this artifact as run input.
            **kwargs: Keyword arguments for the accessor, i.e. `h5py` or `zarr` connection,
                `pyarrow.dataset.dataset`, `polars.scan_*` function.

        Notes:
            For more info, see guide: :doc:`/arrays`.

        Example::

            import lamindb as ln

            # Read AnnData in backed mode from cloud

            artifact = ln.Artifact.get(key="lndb-storage/pbmc68k.h5ad")
            artifact.open()
            #> AnnDataAccessor object with n_obs × n_vars = 70 × 765
            #>     constructed for the AnnData object pbmc68k.h5ad
            #>     ...
            artifact = ln.Artifact.get(key="lndb-storage/df.parquet")
            artifact.open()
            #> pyarrow._dataset.FileSystemDataset

        """
        if self._overwrite_versions and not self.is_latest:
            raise ValueError(INCONSISTENT_STATE_MSG)
        # all hdf5 suffixes including gzipped
        h5_suffixes = [".h5", ".hdf5", ".h5ad"]
        h5_suffixes += [s + ".gz" for s in h5_suffixes]
        # ignore empty suffix for now
        df_suffixes = tuple(set(PYARROW_SUFFIXES).union(POLARS_SUFFIXES))
        suffixes = (
            (
                "",
                ".zarr",
                ".anndata.zarr",
                ".tiledbsoma",
            )
            + tuple(h5_suffixes)
            + df_suffixes
            + tuple(
                s + ".gz" for s in PYARROW_SUFFIXES
            )  # this doesn't work for externally gzipped files, REMOVE LATER
        )
        suffix = self.suffix
        if suffix not in suffixes:
            raise ValueError(
                "Artifact should have a zarr, h5, tiledbsoma object"
                " or a compatible `pyarrow.dataset.dataset` or `polars.scan_*` directory"
                " as the underlying data, please use one of the following suffixes"
                f" for the object name: {', '.join(suffixes[1:])}."
                f" Or no suffix for a folder with {', '.join(df_suffixes)} files"
                " (no mixing allowed)."
            )
        using_key = settings._using_key
        filepath, cache_key = filepath_cache_key_from_artifact(
            self, using_key=using_key
        )

        is_tiledbsoma_w = (
            filepath.name == "soma" or suffix == ".tiledbsoma"
        ) and mode == "w"
        is_zarr_w = suffix == ".zarr" and mode == "r+"

        if mode != "r" and not (is_tiledbsoma_w or is_zarr_w):
            raise ValueError(
                f"It is not allowed to open a {suffix} object with mode='{mode}'. "
                "You can open all supported formats with mode='r', "
                "a tiledbsoma store with mode='w', "
                "AnnData or SpatialData zarr store with mode='r+'."
            )
        # consider the case where an object is already locally cached
        localpath = setup_settings.paths.cloud_to_local_no_update(
            filepath, cache_key=cache_key
        )
        if is_tiledbsoma_w or is_zarr_w:
            open_cache = False
        else:
            open_cache = not isinstance(
                filepath, LocalPathClasses
            ) and not filepath.synchronize_to(localpath, just_check=True)
        if open_cache:
            try:
                access = backed_access(
                    localpath, mode, engine, using_key=using_key, **kwargs
                )
            except Exception as e:
                # also ignore ValueError here because
                # such errors most probably just imply an incorrect argument
                if isinstance(e, (ImportError, ValueError)) or isinstance(
                    filepath, LocalPathClasses
                ):
                    raise e
                logger.warning(
                    f"The cache might be corrupted: {e}. Trying to open directly."
                )
                access = backed_access(
                    filepath, mode, engine, using_key=using_key, **kwargs
                )
                # happens only if backed_access has been successful
                # delete the corrupted cache
                if localpath.is_dir():
                    shutil.rmtree(localpath)
                else:
                    localpath.unlink(missing_ok=True)
        else:
            access = backed_access(self, mode, engine, using_key=using_key, **kwargs)
            if is_tiledbsoma_w:

                def finalize():
                    nonlocal self, filepath, localpath
                    if not isinstance(filepath, LocalPathClasses):
                        _, hash, _, _ = get_stat_dir_cloud(filepath)
                    else:
                        # this can be very slow
                        _, hash, _, _ = hash_dir(filepath)
                    if self.hash != hash:
                        from .sqlrecord import init_self_from_db

                        new_version = Artifact(
                            filepath, revises=self, _is_internal_call=True
                        ).save()
                        # note: sets _state.db = "default"
                        init_self_from_db(self, new_version)

                        if localpath != filepath and localpath.exists():
                            shutil.rmtree(localpath)

                access = _track_writes_factory(access, finalize)
        # only call if open is successfull
        track_run_input(self, is_run_input)
        return access

    def load(
        self, *, is_run_input: bool | None = None, mute: bool = False, **kwargs
    ) -> Any:
        """Cache and load into memory.

        See all :mod:`~lamindb.core.loaders`.

        Args:
            is_run_input: Whether to track this artifact as run input.
            mute: Silence logging of caching progress.
            **kwargs: Keyword arguments for the loader.

        Examples:

            Load a `DataFrame`-like artifact:

            >>> artifact.load().head()
            sepal_length sepal_width petal_length petal_width iris_organism_code
            0        0.051       0.035        0.014       0.002                 0
            1        0.049       0.030        0.014       0.002                 0
            2        0.047       0.032        0.013       0.002                 0
            3        0.046       0.031        0.015       0.002                 0
            4        0.050       0.036        0.014       0.002                 0

            Load an `AnnData`-like artifact:

            >>> artifact.load()
            AnnData object with n_obs × n_vars = 70 × 765

            Fall back to :meth:`~lamindb.Artifact.cache` if no in-memory representation is configured:

            >>> artifact.load()
            PosixPath('/home/runner/work/lamindb/lamindb/docs/guide/mydata/.lamindb/jb7BY5UJoQVGMUOKiLcn.jpg')
        """
        if self._overwrite_versions and not self.is_latest:
            raise ValueError(INCONSISTENT_STATE_MSG)

        if hasattr(self, "_memory_rep") and self._memory_rep is not None:
            access_memory = self._memory_rep
            # SpatialData objects zarr stores are moved when saved
            # SpatialData's __repr__ method attempts to access information from the old path
            # Therefore, we need to update the in-memory path to the now moved Artifact storage path
            if access_memory.__class__.__name__ == "SpatialData":
                access_memory.path = self._cache_path
        else:
            filepath, cache_key = filepath_cache_key_from_artifact(
                self, using_key=settings._using_key
            )
            cache_path = _synchronize_cleanup_on_error(
                filepath, cache_key=cache_key, print_progress=not mute
            )
            try:
                # cache_path is local so doesn't trigger any sync in load_to_memory
                access_memory = load_to_memory(cache_path, **kwargs)
            except Exception as e:
                # raise the exception if it comes from not having a correct loader
                # import error is also most probbaly not a problem with the cache
                # or if the original path is local
                if isinstance(e, (NotImplementedError, ImportError)) or isinstance(
                    filepath, LocalPathClasses
                ):
                    raise e
                logger.warning(
                    f"The cache might be corrupted: {e}. Retrying to synchronize."
                )
                # delete the existing cache
                if cache_path.is_dir():
                    shutil.rmtree(cache_path)
                else:
                    cache_path.unlink(missing_ok=True)
                # download again and try to load into memory
                cache_path = _synchronize_cleanup_on_error(
                    filepath, cache_key=cache_key, print_progress=not mute
                )
                access_memory = load_to_memory(cache_path, **kwargs)
        # only call if load is successfull
        track_run_input(self, is_run_input)

        return access_memory

    def cache(
        self, *, is_run_input: bool | None = None, mute: bool = False, **kwargs
    ) -> Path:
        """Download cloud artifact to local cache.

        Follows synching logic: only caches an artifact if it's outdated in the local cache.

        Returns a path to a locally cached on-disk object (say a `.jpg` file).

        Args:
            mute: Silence logging of caching progress.
            is_run_input: Whether to track this artifact as run input.

        Example::

            # Sync file from cloud and return the local path of the cache
            artifact.cache()
            #> PosixPath('/home/runner/work/Caches/lamindb/lamindb-ci/lndb-storage/pbmc68k.h5ad')
        """
        if self._overwrite_versions and not self.is_latest:
            raise ValueError(INCONSISTENT_STATE_MSG)

        filepath, cache_key = filepath_cache_key_from_artifact(
            self, using_key=settings._using_key
        )
        if mute:
            kwargs["print_progress"] = False
        cache_path = _synchronize_cleanup_on_error(
            filepath, cache_key=cache_key, **kwargs
        )
        # only call if sync is successfull
        track_run_input(self, is_run_input)
        return cache_path

    def delete(
        self,
        permanent: bool | None = None,
        storage: bool | None = None,
        using_key: str | None = None,
    ) -> None:
        """Trash or permanently delete.

        A first call to `.delete()` puts an artifact into the trash (sets `branch_id` to `-1`).
        A second call permanently deletes the artifact.
        If it is a folder artifact with multiple versions, deleting a non-latest version
        will not delete the underlying storage by default (if `storage=True` is not specified).
        Deleting the latest version will delete all the versions for folder artifacts.

        Args:
            permanent: Permanently delete the artifact (skip trash).
            storage: Indicate whether you want to delete the artifact in storage.

        Example::

            import lamindb as ln

            # For an `Artifact` object `artifact`, call:
            artifact = ln.Artifact.get(key="some.csv")
            artifact.delete() # delete a single file artifact

            artifact = ln.Artifact.filter(key="some.tiledbsoma". is_latest=False).first()
            artiact.delete() # delete an old version, the data will not be deleted

            artifact = ln.Artifact.get(key="some.tiledbsoma". is_latest=True)
            artiact.delete() # delete all versions, the data will be deleted or prompted for deletion.
        """
        super().delete(permanent=permanent, storage=storage, using_key=using_key)

    @property
    def _is_saved_to_storage_location(self) -> bool | None:
        """Indicates whether this artifact was correctly written to its storage.

        This is meaningful only after calling `.save()`.

        `None` means no writing was necessary, `True` - that it was written correctly.
        `False` shows that there was a problem with writing.
        """
        if self._aux is not None:
            return self._aux.get("af", {}).get("0", None)
        else:
            return None

    @_is_saved_to_storage_location.setter
    def _is_saved_to_storage_location(self, value: bool) -> None:
        self._aux = self._aux or {}
        self._aux.setdefault("af", {})["0"] = value

    def save(
        self,
        upload: bool | None = None,
        transfer: Literal["record", "annotations"] = "record",
        **kwargs,
    ) -> Artifact:
        """Save to database & storage.

        Args:
            upload: Trigger upload to cloud storage in instances with hybrid storage mode.
            transfer: In case artifact was queried on a different instance, dictates behavior of transfer.
                If "record", only the artifact record is transferred to the current instance.
                If "annotations", also the annotations linked in the source instance are transferred.

        See Also:
            :doc:`transfer`

        Example:

            ::

                import lamindb as ln

                artifact = ln.Artifact("./myfile.csv", key="myfile.parquet").save()
        """
        if transfer not in {"record", "annotations"}:
            raise ValueError(
                f"transfer should be either 'record' or 'annotations', not {transfer}"
            )
        else:
            kwargs["transfer"] = transfer
        state_was_adding = self._state.adding
        print_progress = kwargs.pop("print_progress", True)
        store_kwargs = kwargs.pop(
            "store_kwargs", {}
        )  # kwargs for .upload_from in the end
        access_token = kwargs.pop("access_token", None)
        local_path = None
        if upload and setup_settings.instance.keep_artifacts_local:
            # switch local storage location to cloud
            local_path = self.path
            self.storage_id = setup_settings.instance.storage._id
            self._local_filepath = local_path
            # switch to virtual storage key upon upload
            # the local filepath is already cached at that point
            self._key_is_virtual = True
            # ensure that the artifact is uploaded
            self._to_store = True

        # _is_saved_to_storage_location indicates whether the saving / upload process is successful
        flag_complete = hasattr(self, "_local_filepath") and getattr(
            self, "_to_store", False
        )
        if flag_complete:
            self._is_saved_to_storage_location = (
                False  # will be updated to True at the end
            )

        self._save_skip_storage(**kwargs)

        using_key = None
        if "using" in kwargs:
            using_key = kwargs["using"]
        exception_upload = check_and_attempt_upload(
            self,
            using_key,
            access_token=access_token,
            print_progress=print_progress,
            **store_kwargs,
        )
        if exception_upload is not None:
            # we do not want to raise file not found on cleanup if upload of a file failed
            # often it is ACID in the filesystem itself
            # for example, s3 won't have the failed file, so just skip the delete in this case
            raise_file_not_found_error = False
            self._delete_skip_storage()
        else:
            # this is the case when it is cleaned on .replace
            raise_file_not_found_error = True
        # this is triggered by an exception in check_and_attempt_upload or by replace.
        exception_clear = check_and_attempt_clearing(
            self,
            raise_file_not_found_error=raise_file_not_found_error,
            using_key=using_key,
        )
        if exception_upload is not None:
            raise exception_upload
        if exception_clear is not None:
            raise exception_clear
        # the saving / upload process has been successfull, just mark it as such
        # maybe some error handling here?
        if flag_complete:
            self._is_saved_to_storage_location = True
            # pass kwargs here because it can contain `using` or other things
            # affecting the connection
            super().save(**kwargs)

        # this is only for keep_artifacts_local
        if local_path is not None and not state_was_adding:
            # only move the local artifact to cache if it was not newly created
            local_path_cache = ln_setup.settings.cache_dir / local_path.name
            # don't use Path.rename here because of cross-device link error
            # https://laminlabs.slack.com/archives/C04A0RMA0SC/p1710259102686969
            shutil.move(
                local_path,  # type: ignore
                local_path_cache,
            )
            logger.important(f"moved local artifact to cache: {local_path_cache}")

        # annotate with external features
        if hasattr(self, "_external_features"):
            external_features = self._external_features
            self.features.set_values(external_features)
        # annotate with internal features based on curator
        if hasattr(self, "_curator"):
            curator = self._curator
            delattr(self, "_curator")
            # just annotates this artifact
            curator.save_artifact()
        if hasattr(self, "_external_features"):
            delattr(self, "_external_features")
        return self

    def describe(self, return_str: bool = False) -> None | str:
        """Describe record including relations.

        Args:
            return_str: Return a string instead of printing.
        """
        from ._describe import describe_postgres_sqlite

        return describe_postgres_sqlite(self, return_str=return_str)


# can't really just call .cache in .load because of double tracking
def _synchronize_cleanup_on_error(
    filepath: UPath, cache_key: str | None = None, **kwargs
) -> UPath:
    try:
        print_progress = kwargs.pop("print_progress", True)
        cache_path = setup_settings.paths.cloud_to_local(
            filepath, cache_key=cache_key, print_progress=print_progress, **kwargs
        )
    except Exception as e:
        if not isinstance(filepath, LocalPathClasses):
            cache_path = setup_settings.paths.cloud_to_local_no_update(
                filepath, cache_key=cache_key
            )
            if cache_path.is_dir():
                shutil.rmtree(cache_path)
            else:
                cache_path.unlink(missing_ok=True)
        raise e
    return cache_path


def _delete_skip_storage(artifact, *args, **kwargs) -> None:
    super(SQLRecord, artifact).delete(*args, **kwargs)


def _save_skip_storage(artifact, **kwargs) -> None:
    save_staged_feature_sets(artifact)
    super(Artifact, artifact).save(**kwargs)
    save_schema_links(artifact)


class ArtifactFeatureValue(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(
        Artifact, CASCADE, related_name="links_featurevalue"
    )
    # we follow the lower() case convention rather than snake case for link models
    featurevalue = ForeignKey(FeatureValue, PROTECT, related_name="links_artifact")

    class Meta:
        app_label = "lamindb"
        unique_together = ("artifact", "featurevalue")


class ArtifactUser(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey("Artifact", CASCADE, related_name="links_user")
    user: User = ForeignKey(User, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_artifactuser", default=None
    )
    label_ref_is_name: bool | None = BooleanField(null=True)
    feature_ref_is_name: bool | None = BooleanField(null=True)

    class Meta:
        # can have the same label linked to the same artifact if the feature is
        # different
        app_label = "lamindb"
        unique_together = ("artifact", "user", "feature")


def track_run_input(
    record: (
        Artifact | Iterable[Artifact]
    ),  # can also be Collection | Iterable[Collection]
    is_run_input: bool | Run | None = None,
    run: Run | None = None,
) -> None:
    """Links a record as an input to a run.

    This function contains all validation logic to make decisions on whether a
    record qualifies as an input or not.
    """
    if is_run_input is False:
        return None

    from .._tracked import get_current_tracked_run
    from ..core._context import context
    from .collection import Collection

    if isinstance(is_run_input, Run):
        run = is_run_input
        is_run_input = True
    elif run is None:
        run = get_current_tracked_run()
        if run is None:
            run = context.run
    # consider that record is an iterable of Data
    record_iter: Iterable[Artifact] | Iterable[Collection] = (
        [record] if isinstance(record, (Artifact, Collection)) else record
    )
    input_records = []
    if run is not None:
        assert not run._state.adding, "Save the run before tracking its inputs."  # noqa: S101

        def is_valid_input(record: Artifact | Collection):
            is_valid = False
            # if a record is not yet saved it has record._state.db = None
            # then it can't be an input
            # we silently ignore because what will happen is that
            # the record either gets saved and then is tracked as an output
            # or it won't get saved at all
            if record._state.db == "default":
                # things are OK if the record is on the default db
                is_valid = True
            else:
                # record is on another db
                # we have to save the record into the current db with
                # the run being attached to a transfer transform
                logger.info(
                    f"completing transfer to track {record.__class__.__name__}('{record.uid}') as input"
                )
                record.save()
                is_valid = True
            # avoid cycles: record can't be both input and output
            if record.run_id == run.id:
                logger.debug(
                    f"not tracking {record} as input to run {run} because created by same run"
                )
                is_valid = False
            if run.id == getattr(record, "_subsequent_run_id", None):
                logger.debug(
                    f"not tracking {record} as input to run {run} because re-created in same run"
                )
                is_valid = False
            return is_valid

        input_records = [record for record in record_iter if is_valid_input(record)]
        input_records_ids = [record.id for record in input_records]
    if input_records:
        record_class_name = input_records[0].__class__.__name__.lower()
    # let us first look at the case in which the user does not
    # provide a boolean value for `is_run_input`
    # hence, we need to determine whether we actually want to
    # track a run or not
    track = False
    is_run_input = settings.track_run_inputs if is_run_input is None else is_run_input
    if is_run_input:
        if run is None:
            if not is_read_only_connection():
                logger.warning(WARNING_NO_INPUT)
        elif input_records:
            logger.debug(
                f"adding {record_class_name} ids {input_records_ids} as inputs for run {run.id}"
            )
            track = True
    else:
        track = is_run_input
    if not track or not input_records:
        return None
    if run is None:
        raise ValueError("No run context set. Call `ln.track()`.")
    if record_class_name == "artifact":
        IsLink = run.input_artifacts.through
        links = [
            IsLink(run_id=run.id, artifact_id=record_id)
            for record_id in input_records_ids
        ]
    else:
        IsLink = run.input_collections.through
        links = [
            IsLink(run_id=run.id, collection_id=record_id)
            for record_id in input_records_ids
        ]
    try:
        IsLink.objects.bulk_create(links, ignore_conflicts=True)
    except ProgrammingError as e:
        if "new row violates row-level security policy" in str(e):
            instance = setup_settings.instance
            available_spaces = instance.available_spaces
            if available_spaces is None:
                raise NoWriteAccess(
                    f"You’re not allowed to write to the instance {instance.slug}.\n"
                    "Please contact administrators of the instance if you need write access."
                ) from None
            write_access_spaces = available_spaces["admin"] + available_spaces["write"]
            no_write_access_spaces = {
                record_space
                for record in input_records
                if (record_space := record.space) not in write_access_spaces
            }
            if (run_space := run.space) not in write_access_spaces:
                no_write_access_spaces.add(run_space)

            if not no_write_access_spaces:
                # if there are no unavailable spaces, then this should be due to locking
                locked_records = [
                    record
                    for record in input_records
                    if getattr(record, "is_locked", False)
                ]
                if run.is_locked:
                    locked_records.append(run)
                # if no unavailable spaces and no locked records, just raise the original error
                if not locked_records:
                    raise e
                no_write_msg = (
                    "It is not allowed to modify locked records: "
                    + ", ".join(
                        r.__class__.__name__ + f"(uid={r.uid})" for r in locked_records
                    )
                    + "."
                )
                raise NoWriteAccess(no_write_msg) from None

            if len(no_write_access_spaces) > 1:
                name_msg = ", ".join(
                    f"'{space.name}'" for space in no_write_access_spaces
                )
                space_msg = "spaces"
            else:
                name_msg = f"'{no_write_access_spaces.pop().name}'"
                space_msg = "space"
            raise NoWriteAccess(
                f"You’re not allowed to write to the {space_msg} {name_msg}.\n"
                f"Please contact administrators of the {space_msg} if you need write access."
            ) from None
        else:
            raise e


# privates currently dealt with separately
# mypy: ignore-errors
Artifact._delete_skip_storage = _delete_skip_storage
Artifact._save_skip_storage = _save_skip_storage
Artifact.view_lineage = view_lineage
