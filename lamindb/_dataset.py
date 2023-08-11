from typing import Iterable, List, Optional, Union

import anndata as ad
import pandas as pd
from lnschema_core import ids
from lnschema_core.models import Dataset

from . import FeatureSet, File, Run
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
    data: Optional[Union[pd.DataFrame, ad.AnnData]] = None
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
        if isinstance(data, pd.DataFrame):
            # description filled in below
            file = File.from_df(data, run=run, description="See dataset {}")
            dataset._feature_sets = file._feature_sets
        elif isinstance(data, ad.AnnData):
            if len(feature_sets) != 2:
                raise ValueError(
                    "Please provide a feature set describing each `.var.index` &"
                    " `.obs.columns`"
                )
            file = File.from_anndata(
                data, run=run, feature_sets=feature_sets, description="See dataset {}"
            )
            dataset._feature_sets = feature_sets
        else:
            raise ValueError("Only DataFrame and AnnData can be passed as data")
        hash = file.hash
        id = file.id
        file.description = f"See dataset {id}"
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
    # query all feature sets at the same time rather than making a single query per file
    feature_set_file_links = File.feature_sets.through.objects.filter(
        file_id__in=file_ids
    )
    feature_set_ids = [link.feature_set_id for link in feature_set_file_links]
    feature_sets = FeatureSet.filter(id__in=feature_set_ids)
    # validate consistency of hashes
    # we do not allow duplicate hashes
    file_hashes = [file.hash for file in files]
    file_hashes_set = set(file_hashes)
    if len(file_hashes) != len(file_hashes_set):
        raise ValueError("Please pass distinct files")
    hash = hash_set(file_hashes_set)
    # create the dataset
    dataset = Dataset(name=name, hash=hash, feature_sets=feature_sets, files=files)
    return dataset


def backed(dataset: Dataset):
    if dataset.file is None:
        raise RuntimeError("Can only call backed() for datasets with a single file")
    return dataset.file.backed()


def load(dataset: Dataset):
    """Load the combined dataset."""
    if dataset.file is not None:
        return dataset.file.load()
    else:
        suffixes = [file.suffix for file in dataset.files.all()]
        if len(set(suffixes)) != 1:
            raise RuntimeError(
                "Can only load datasets where all files have the same suffix"
            )
        objects = [file.load() for file in dataset.files.all()]
        if isinstance(objects[0], pd.DataFrame):
            return pd.concat(objects)
        elif isinstance(objects[0], ad.AnnData):
            return ad.concat(objects)


def delete(dataset: Dataset, storage: bool = False):
    super(Dataset, dataset).delete()
    if dataset.file is not None:
        dataset.file.delete(storage=storage)


def save(dataset: Dataset):
    if dataset.file is not None:
        dataset.file.save()
    feature_sets = dataset._feature_sets
    if isinstance(feature_sets, dict):
        feature_sets = feature_sets.values()
    for feature_set in feature_sets:
        feature_set.save()
    super(Dataset, dataset).save()
    if len(dataset._files) > 0:
        dataset.files.set(dataset._files)
    if len(dataset._feature_sets) > 0:
        dataset.feature_sets.set(feature_sets)


Dataset.__init__ = __init__
Dataset.from_files = from_files
Dataset.backed = backed
Dataset.load = load
Dataset.delete = delete
Dataset.save = save
