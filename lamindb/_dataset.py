from typing import Iterable, List, Optional, Union

import pandas as pd
from anndata import AnnData
from lnschema_core import ids
from lnschema_core.models import Dataset
from pandas import DataFrame

from . import Feature, FeatureSet, File, Run
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
    data: Optional[Union[DataFrame, AnnData]] = None
    if "data" in kwargs or len(args) == 1:
        data = kwargs.pop("data") if len(args) == 0 else args[0]
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    run: Optional[Run] = kwargs.pop("run") if "run" in kwargs else None
    files: List[File] = kwargs.pop("files") if "files" in kwargs else []
    file: Optional[File] = kwargs.pop("file") if "file" in kwargs else None
    hash: Optional[str] = kwargs.pop("hash") if "hash" in kwargs else None
    feature_sets: List[FeatureSet] = (
        kwargs.pop("feature_sets") if "feature_sets" in kwargs else []
    )
    assert len(kwargs) == 0
    if data is not None:
        if isinstance(data, DataFrame):
            feature_set = FeatureSet.from_values(data.columns, Feature.name)
        file = File(data=data, name=name, run=run)
        dataset._feature_sets = [feature_set]
        hash = file.hash
        id = file.id
    else:
        id = ids.base62_20()
        dataset._feature_sets = feature_sets
    super(Dataset, dataset).__init__(id=id, name=name, file=file, hash=hash)
    dataset._files = files


@classmethod  # type: ignore
def from_files(dataset: Dataset, *, name: str, files: Iterable[File]) -> Dataset:
    # assert all files are already saved
    # saved = not any([file._state._adding for file in files])
    # if not saved:
    #     raise ValueError("Not all files are yet saved, please save them")
    # query all feature sets of files
    file_ids = [file.id for file in files]
    feature_sets = list(
        File.feature_sets.through.objects.filter(file_id__in=file_ids).all()
    )
    # validate consistency of feature_sets
    # we only allow one feature set per type
    feature_set_types = [feature_set.type for feature_set in feature_sets]
    feature_set_ids_types = [
        (feature_set.id, feature_set.type) for feature_set in feature_sets
    ]
    if len(set(feature_set_ids_types)) != len(set(feature_set_types)):
        # we can do below in the future!
        # logger.warning(
        #     "feature sets are inconsistent across files"
        #     "computing union! files will be outer-joined"
        # )
        raise ValueError(
            "Currently only supporting datasets from files with same feature sets"
        )
    # validate consistency of hashes
    # we do not allow duplicate hashes
    file_hashes = [file.hash for file in files]
    file_hashes_set = set(file_hashes)
    assert len(file_hashes) == len(file_hashes_set)
    hash = hash_set(file_hashes_set)
    dataset = Dataset(name=name, hash=hash, feature_sets=feature_sets, files=files)
    return dataset


def backed(dataset: Dataset):
    return dataset.file.backed()


def load(dataset: Dataset):
    """Load the combined dataset."""
    if dataset.file is not None:
        return dataset.file.load()
    else:
        objects = [file.load() for file in dataset.files.all()]
        if isinstance(objects[0], pd.DataFrame):
            return pd.concat(objects)


def delete(dataset: Dataset):
    dataset.file.delete()
    super(Dataset, dataset).delete()


def save(dataset: Dataset):
    if dataset.file is not None:
        dataset.file.save()
    for feature_set in dataset._feature_sets:
        feature_set.save()
    super(Dataset, dataset).save()
    if len(dataset._files) > 0:
        dataset.files.set(dataset._files)
    if len(dataset._feature_sets) > 0:
        dataset.feature_sets.set(dataset._feature_sets)


Dataset.__init__ = __init__
Dataset.from_files = from_files
Dataset.backed = backed
Dataset.load = load
Dataset.delete = delete
Dataset.save = save
