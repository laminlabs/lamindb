from typing import Dict, Iterable, Optional, Tuple, Union

import anndata as ad
import pandas as pd
from lnschema_core.models import Dataset

from . import File, Run
from .dev._data import add_transform_to_kwargs, get_run, save_feature_set_links
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
    description: Optional[str] = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    run: Optional[Run] = kwargs.pop("run") if "run" in kwargs else None
    assert len(kwargs) == 0
    id = None
    feature_sets = None
    run = get_run(run)
    # init file
    if isinstance(data, pd.DataFrame) or isinstance(data, ad.AnnData):
        files = None
        if isinstance(data, pd.DataFrame):
            file = File.from_df(data, run=run, description="tmp")
        elif isinstance(data, ad.AnnData):
            file = File.from_anndata(data, run=run, description="tmp")
        feature_sets = file._feature_sets  # type: ignore
        hash = file.hash  # type: ignore
        id = file.id  # type: ignore
        file.description = f"See dataset {id}"  # type: ignore
    # init files
    else:
        file = None
        if hasattr(data, "__getitem__"):
            assert isinstance(data[0], File)  # type: ignore
            files = data
            hash, feature_sets = from_files(files)  # type: ignore
        else:
            raise ValueError("Only DataFrame, AnnData and iterable of File is allowed")
    kwargs = {}
    add_transform_to_kwargs(kwargs, run)
    super(Dataset, dataset).__init__(
        id=id,
        name=name,
        description=description,
        file=file,
        hash=hash,
        run=run,
        **kwargs,
    )
    dataset._files = files
    dataset._feature_sets = feature_sets


def from_files(files: Iterable[File]) -> Tuple[str, Dict[str, str]]:
    # assert all files are already saved
    saved = not any([file._state.adding for file in files])
    if not saved:
        raise ValueError("Not all files are yet saved, please save them")
    # query all feature sets of files
    file_ids = [file.id for file in files]
    # query all feature sets at the same time rather than making a single query per file
    feature_set_file_links = File.feature_sets.through.objects.filter(
        file_id__in=file_ids
    )
    feature_set_slots_ids = {}
    for link in feature_set_file_links:
        feature_set_slots_ids[link.slot] = link.feature_set_id
    # validate consistency of hashes
    # we do not allow duplicate hashes
    hashes = [file.hash for file in files]
    if len(hashes) != len(set(hashes)):
        seen = set()
        non_unique = [x for x in hashes if x in seen or seen.add(x)]  # type: ignore
        raise ValueError(
            "Please pass files with distinct hashes: these ones are non-unique"
            f" {non_unique}"
        )
    hash = hash_set(set(hashes))
    return hash, feature_set_slots_ids


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
    # we don't need to save feature sets again
    # save_transform_run_feature_sets(dataset)
    super(Dataset, dataset).save()
    if dataset._files is not None and len(dataset._files) > 0:
        dataset.files.set(dataset._files)
    save_feature_set_links(dataset)


Dataset.__init__ = __init__
Dataset.from_files = from_files
Dataset.backed = backed
Dataset.load = load
Dataset.delete = delete
Dataset.save = save
