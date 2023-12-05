from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Optional, Tuple, Union

import anndata as ad
import pandas as pd
from lamin_utils import logger
from lamindb_setup._init_instance import register_storage
from lamindb_setup.dev import StorageSettings
from lamindb_setup.dev._docs import doc_args
from lamindb_setup.dev._hub_utils import get_storage_region
from lamindb_setup.dev.upath import UPath
from lnschema_core.models import Dataset, Feature, FeatureSet
from lnschema_core.types import AnnDataLike, DataLike, FieldAttr, VisibilityChoice

from lamindb._utils import attach_func_to_class_method
from lamindb.dev._data import _track_run_input
from lamindb.dev._mapped_dataset import MappedDataset
from lamindb.dev.storage._backed_access import AnnDataAccessor, BackedAccessor
from lamindb.dev.versioning import get_ids_from_old_version, init_uid

from . import _TESTING, File, Run
from ._file import parse_feature_sets_from_anndata
from ._registry import init_self_from_db
from .dev._data import (
    add_transform_to_kwargs,
    get_run,
    save_feature_set_links,
    save_feature_sets,
)
from .dev.hashing import hash_set


def __init__(
    dataset: Dataset,
    *args,
    **kwargs,
):
    if len(args) == len(dataset._meta.concrete_fields):
        super(Dataset, dataset).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: data")
    data: Union[pd.DataFrame, ad.AnnData, File, Iterable[File]] = (
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
    is_new_version_of: Optional[Dataset] = (
        kwargs.pop("is_new_version_of") if "is_new_version_of" in kwargs else None
    )
    initial_version_id: Optional[int] = (
        kwargs.pop("initial_version_id") if "initial_version_id" in kwargs else None
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
    if not len(kwargs) == 0:
        raise ValueError(
            f"Only data, name, run, description, reference, reference_type, visibility can be passed, you passed: {kwargs}"  # noqa
        )

    if is_new_version_of is None:
        provisional_uid = init_uid(version=version, n_full_id=20)
    else:
        if not isinstance(is_new_version_of, Dataset):
            raise TypeError("is_new_version_of has to be of type ln.Dataset")
        provisional_uid, initial_version_id, version = get_ids_from_old_version(
            is_new_version_of, version, n_full_id=20
        )
        if name is None:
            name = is_new_version_of.name
    if version is not None:
        if initial_version_id is None:
            logger.info(
                "initializing versioning for this dataset! create future versions of it"
                " using ln.Dataset(..., is_new_version_of=old_dataset)"
            )

    run = get_run(run)
    data_init_complete = False
    file = None
    files = None
    storage = None
    # init from directory or bucket
    if isinstance(data, (str, Path, UPath)):
        upath = UPath(data)
        # below frequently times out on GCP
        # comment this and corresponding test out
        # if not upath.is_dir():
        #     raise ValueError(f"Can only pass buckets or directories, not {data}")
        upath_str = upath.as_posix().rstrip("/")
        region = get_storage_region(upath_str)
        storage_settings = StorageSettings(upath_str, region)
        storage = register_storage(storage_settings)
        hash = None
        data_init_complete = True
    # now handle potential metadata
    if meta is not None:
        if not isinstance(meta, (pd.DataFrame, ad.AnnData, File)):
            raise ValueError(
                "meta has to be of type `(pd.DataFrame, ad.AnnData, File)`"
            )
        data = meta
    # init file - is either data or metadata
    if isinstance(data, (pd.DataFrame, ad.AnnData, File)):
        if isinstance(data, File):
            file = data
            if file._state.adding:
                raise ValueError("Save file before creating dataset!")
            if not feature_sets:
                feature_sets = file.features._feature_set_by_slot
            else:
                if len(file.features._feature_set_by_slot) > 0:
                    logger.info("overwriting feature sets linked to file")
        else:
            log_hint = True if feature_sets is None else False
            file_is_new_version_of = (
                is_new_version_of.file if is_new_version_of is not None else None
            )
            file = File(
                data,
                run=run,
                description="tmp",
                log_hint=log_hint,
                version=version,
                is_new_version_of=file_is_new_version_of,
            )
            # do we really want to update the file here?
            if feature_sets:
                file._feature_sets = feature_sets
        hash = file.hash  # type: ignore
        provisional_uid = file.uid  # type: ignore
        if file.description is None or file.description == "tmp":
            file.description = f"See dataset {provisional_uid}"  # type: ignore
        data_init_complete = True
    if not data_init_complete:
        if hasattr(data, "__getitem__"):
            assert isinstance(data[0], File)  # type: ignore
            files = data
            hash, feature_sets = from_files(files)  # type: ignore
            data_init_complete = True
        else:
            raise ValueError(
                "Only DataFrame, AnnData, folder or list of File is allowed."
            )
    # we ignore datasets in trash containing the same hash
    if hash is not None:
        existing_dataset = Dataset.filter(hash=hash).one_or_none()
    else:
        existing_dataset = None
    if existing_dataset is not None:
        logger.warning(f"returning existing dataset with same hash: {existing_dataset}")
        init_self_from_db(dataset, existing_dataset)
        for slot, feature_set in dataset.features._feature_set_by_slot.items():
            if slot in feature_sets:
                if not feature_sets[slot] == feature_set:
                    dataset.feature_sets.remove(feature_set)
                    logger.warning(f"removing feature set: {feature_set}")
    else:
        kwargs = {}
        add_transform_to_kwargs(kwargs, run)
        super(Dataset, dataset).__init__(
            uid=provisional_uid,
            name=name,
            description=description,
            reference=reference,
            reference_type=reference_type,
            file=file,
            storage=storage,
            hash=hash,
            run=run,
            version=version,
            initial_version_id=initial_version_id,
            visibility=visibility,
            **kwargs,
        )
    dataset._files = files
    dataset._feature_sets = feature_sets
    # register provenance
    if is_new_version_of is not None:
        _track_run_input(is_new_version_of, run=run)
    if file is not None and file.run != run:
        _track_run_input(file, run=run)
    elif files is not None:
        _track_run_input(files, run=run)


@classmethod  # type: ignore
@doc_args(Dataset.from_df.__doc__)
def from_df(
    cls,
    df: "pd.DataFrame",
    field: FieldAttr = Feature.name,
    name: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
    reference: Optional[str] = None,
    reference_type: Optional[str] = None,
    version: Optional[str] = None,
    is_new_version_of: Optional["File"] = None,
    **kwargs,
) -> "Dataset":
    """{}"""
    feature_set = FeatureSet.from_df(df, field=field, **kwargs)
    if feature_set is not None:
        feature_sets = {"columns": feature_set}
    else:
        feature_sets = {}
    dataset = Dataset(
        data=df,
        name=name,
        run=run,
        description=description,
        feature_sets=feature_sets,
        reference=reference,
        reference_type=reference_type,
        version=version,
        is_new_version_of=is_new_version_of,
    )
    return dataset


@classmethod  # type: ignore
@doc_args(Dataset.from_anndata.__doc__)
def from_anndata(
    cls,
    adata: "AnnDataLike",
    field: Optional[FieldAttr],
    name: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
    reference: Optional[str] = None,
    reference_type: Optional[str] = None,
    version: Optional[str] = None,
    is_new_version_of: Optional["File"] = None,
    **kwargs,
) -> "Dataset":
    """{}"""
    if isinstance(adata, File):
        assert not adata._state.adding
        assert adata.accessor == "AnnData"
        adata_parse = adata.path
    else:
        adata_parse = adata
    feature_sets = parse_feature_sets_from_anndata(adata_parse, field, **kwargs)
    dataset = Dataset(
        data=adata,
        run=run,
        name=name,
        description=description,
        feature_sets=feature_sets,
        reference=reference,
        reference_type=reference_type,
        version=version,
        is_new_version_of=is_new_version_of,
    )
    return dataset


# internal function, not exposed to user
def from_files(files: Iterable[File]) -> Tuple[str, Dict[str, str]]:
    # assert all files are already saved
    logger.debug("check not saved")
    saved = not any([file._state.adding for file in files])
    if not saved:
        raise ValueError("Not all files are yet saved, please save them")
    # query all feature sets of files
    logger.debug("file ids")
    file_ids = [file.id for file in files]
    # query all feature sets at the same time rather than making a single query per file
    logger.debug("feature_set_file_links")
    feature_set_file_links = File.feature_sets.through.objects.filter(
        file_id__in=file_ids
    )
    feature_sets_by_slots = defaultdict(list)
    logger.debug("slots")
    for link in feature_set_file_links:
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
    # file.hash is None for zarr
    # todo: more careful handling of such cases
    hashes = [file.hash for file in files if file.hash is not None]
    if len(hashes) != len(set(hashes)):
        seen = set()
        non_unique = [x for x in hashes if x in seen or seen.add(x)]  # type: ignore
        raise ValueError(
            "Please pass files with distinct hashes: these ones are non-unique"
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
    join_vars: Optional[Literal["auto", "inner"]] = "auto",
    encode_labels: bool = True,
    parallel: bool = False,
    stream: bool = False,
    is_run_input: Optional[bool] = None,
) -> "MappedDataset":
    _track_run_input(self, is_run_input)
    path_list = []
    for file in self.files.all():
        if file.suffix not in {".h5ad", ".zrad", ".zarr"}:
            logger.warning(f"Ignoring file with suffix {file.suffix}")
            continue
        elif not stream and file.suffix == ".h5ad":
            path_list.append(file.stage())
        else:
            path_list.append(file.path)
    return MappedDataset(path_list, label_keys, join_vars, encode_labels, parallel)


# docstring handled through attach_func_to_class_method
def backed(
    self, is_run_input: Optional[bool] = None
) -> Union["AnnDataAccessor", "BackedAccessor"]:
    _track_run_input(self, is_run_input)
    if self.file is None:
        raise RuntimeError("Can only call backed() for datasets with a single file")
    return self.file.backed()


# docstring handled through attach_func_to_class_method
def load(
    self,
    join: Literal["inner", "outer"] = "outer",
    is_run_input: Optional[bool] = None,
    **kwargs,
) -> DataLike:
    # cannot call _track_run_input here, see comment further down
    if self.file is not None:
        _track_run_input(self, is_run_input)
        return self.file.load()
    else:
        all_files = self.files.all()
        suffixes = [file.suffix for file in all_files]
        if len(set(suffixes)) != 1:
            raise RuntimeError(
                "Can only load datasets where all files have the same suffix"
            )
        # because we're tracking data flow on the dataset-level, here, we don't
        # want to track it on the file-level
        objects = [file.load(is_run_input=False) for file in all_files]
        file_uids = [file.uid for file in all_files]
        if isinstance(objects[0], pd.DataFrame):
            concat_object = pd.concat(objects, join=join)
        elif isinstance(objects[0], ad.AnnData):
            concat_object = ad.concat(
                objects, join=join, label="file_uid", keys=file_uids
            )
        # only call it here because there might be errors during concat
        _track_run_input(self, is_run_input)
        return concat_object


# docstring handled through attach_func_to_class_method
def delete(
    self, permanent: Optional[bool] = None, storage: Optional[bool] = None
) -> None:
    # change visibility to trash
    if self.visibility > VisibilityChoice.trash.value and permanent is not True:
        self.visibility = VisibilityChoice.trash.value
        self.save()
        logger.warning("moved dataset to trash.")
        if self.file is not None:
            self.file.visibility = VisibilityChoice.trash.value
            self.file.save()
            logger.warning("moved dataset.file to trash.")
        return

    # permanent delete
    if permanent is None:
        response = input(
            "Dataset record is already in trash! Are you sure to delete it from your"
            " database? (y/n) You can't undo this action."
        )
        delete_record = response == "y"
    else:
        delete_record = permanent

    if delete_record:
        super(Dataset, self).delete()
    if self.file is not None:
        self.file.delete(permanent=permanent, storage=storage)


# docstring handled through attach_func_to_class_method
def save(self, *args, **kwargs) -> None:
    if self.file is not None:
        self.file.save()
    # we don't need to save feature sets again
    save_feature_sets(self)
    super(Dataset, self).save()
    if hasattr(self, "_files"):
        if self._files is not None and len(self._files) > 0:
            self.files.set(self._files)
    save_feature_set_links(self)


@property  # type: ignore
@doc_args(Dataset.path.__doc__)
def path(self) -> Union[Path, UPath]:
    """{}"""
    _track_run_input(self)
    return self.storage.path


# docstring handled through attach_func_to_class_method
def restore(self) -> None:
    self.visibility = VisibilityChoice.default.value
    self.save()
    if self.file is not None:
        self.file.visibility = VisibilityChoice.default.value
        self.file.save()


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

if _TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(Dataset, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Dataset, globals())

setattr(Dataset, "path", path)
# this seems a Django-generated function
delattr(Dataset, "get_visibility_display")
