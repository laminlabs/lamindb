from typing import Optional, Union

from anndata import AnnData
from lnschema_core.models import Dataset
from pandas import DataFrame

from . import Feature, FeatureSet, File, Run


def __init__(
    dataset: Dataset,
    *args,
    **kwargs,
):
    # Below checks for the Django-internal call in from_db()
    # it'd be better if we could avoid this, but not being able to create a File
    # from data with the default constructor renders the central class of the API
    # essentially useless
    # The danger below is not that a user might pass as many args (12 of it), but rather
    # that at some point the Django API might change; on the other hand, this
    # condition of for calling the constructor based on kwargs should always
    # stay robust
    if len(args) == len(dataset._meta.concrete_fields):
        super(Dataset, dataset).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: data")
    data: Union[DataFrame, AnnData] = kwargs.pop("data") if len(args) == 0 else args[0]
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    run: Optional[Run] = kwargs.pop("run") if "run" in kwargs else None
    assert len(kwargs) == 0
    if isinstance(data, DataFrame):
        feature_set = FeatureSet.from_values(data.columns, Feature.name)
    file = File(data=data, name=name, run=run)
    dataset._feature_sets = [feature_set]
    super(Dataset, dataset).__init__(id=file.id, name=name, file=file)


def backed(dataset: Dataset):
    return dataset.file.backed()


def load(dataset: Dataset):
    return dataset.file.load()


def delete(dataset: Dataset):
    dataset.file.delete()
    super(Dataset, dataset).delete()


def save(dataset: Dataset):
    if dataset.file is not None:
        dataset.file.save()
    for feature_set in dataset._feature_sets:
        feature_set.save()
    super(Dataset, dataset).save()


Dataset.__init__ = __init__
Dataset.backed = backed
Dataset.load = load
Dataset.delete = delete
Dataset.save = save
