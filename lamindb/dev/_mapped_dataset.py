from collections import Counter
from os import PathLike
from typing import List, Optional, Union

import numpy as np
from lamindb_setup.dev.upath import UPath

from .storage._backed_access import ArrayTypes, GroupTypes, StorageType, registry


class MappedDataset:
    """Map-style dataset for use in data loaders.

    This currently only works for collections of `AnnData` objects.

    For an example, see :meth:`~lamindb.Dataset.mapped`.

    .. note::

        A similar data loader exists `here
        <https://github.com/Genentech/scimilarity>`__.
    """

    def __init__(
        self,
        path_list: List[Union[str, PathLike]],
        label_keys: Optional[Union[str, List[str]]] = None,
        encode_labels: bool = True,
    ):
        self.storages = []
        self.conns = []
        for path in path_list:
            path = UPath(path)
            if path.exists() and path.is_file():  # type: ignore
                conn, storage = registry.open("h5py", path)
            else:
                conn, storage = registry.open("zarr", path)
            self.conns.append(conn)
            self.storages.append(storage)

        self.n_obs_list = []
        for storage in self.storages:
            X = storage["X"]
            if isinstance(X, ArrayTypes):  # type: ignore
                self.n_obs_list.append(X.shape[0])
            else:
                self.n_obs_list.append(X.attrs["shape"][0])
        self.n_obs = sum(self.n_obs_list)

        self.indices = np.hstack([np.arange(n_obs) for n_obs in self.n_obs_list])
        self.storage_idx = np.repeat(np.arange(len(self.storages)), self.n_obs_list)

        self.encode_labels = encode_labels
        if isinstance(label_keys, str):
            label_keys = [label_keys]
        self.label_keys = label_keys
        if self.label_keys is not None:
            if self.encode_labels:
                self.encoders = []
                for label in self.label_keys:
                    cats = self.get_merged_categories(label)
                    self.encoders.append({cat: i for i, cat in enumerate(cats)})

    def __len__(self):
        return self.n_obs

    def __getitem__(self, idx):
        obs_idx = self.indices[idx]
        storage = self.storages[self.storage_idx[idx]]
        out = [self.get_data_idx(storage, obs_idx)]
        if self.label_keys is not None:
            for i, label in enumerate(self.label_keys):
                label_idx = self.get_label_idx(storage, obs_idx, label)
                if self.encode_labels:
                    label_idx = self.encoders[i][label_idx]
                out.append(label_idx)
        return out

    def get_data_idx(
        self, storage: StorageType, idx: int, layer_key: Optional[str] = None  # type: ignore # noqa
    ):
        """Get the index for the data."""
        layer = storage["X"] if layer_key is None else storage["layers"][layer_key]  # type: ignore # noqa
        if isinstance(layer, ArrayTypes):  # type: ignore
            return layer[idx]
        else:  # assume csr_matrix here
            data = layer["data"]
            indices = layer["indices"]
            indptr = layer["indptr"]
            s = slice(*(indptr[idx : idx + 2]))
            layer_idx = np.zeros(layer.attrs["shape"][1])
            layer_idx[indices[s]] = data[s]
            return layer_idx

    def get_label_idx(self, storage: StorageType, idx: int, label_key: str):  # type: ignore # noqa
        """Get the index for the label by key."""
        obs = storage["obs"]  # type: ignore
        # how backwards compatible do we want to be here actually?
        if isinstance(obs, ArrayTypes):  # type: ignore
            label = obs[idx][obs.dtype.names.index(label_key)]
        else:
            labels = obs[label_key]
            if isinstance(labels, ArrayTypes):  # type: ignore
                label = labels[idx]
            else:
                label = labels["codes"][idx]

        cats = self.get_categories(storage, label_key)
        if cats is not None:
            label = cats[label]
        if isinstance(label, bytes):
            label = label.decode("utf-8")
        return label

    def get_label_weights(self, label_key: str):
        """Get all weights for a given label key."""
        labels = self.get_merged_labels(label_key)
        counter = Counter(labels)  # type: ignore
        weights = np.array([counter[label] for label in labels]) / len(labels)
        return weights

    def get_merged_labels(self, label_key: str):
        """Get merged labels."""
        labels_merge = []
        decode = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)
        for storage in self.storages:
            codes = self.get_codes(storage, label_key)
            labels = decode(codes) if isinstance(codes[0], bytes) else codes
            cats = self.get_categories(storage, label_key)
            if cats is not None:
                cats = decode(cats) if isinstance(cats[0], bytes) else cats
                labels = cats[labels]
            labels_merge.append(labels)
        return np.hstack(labels_merge)

    def get_merged_categories(self, label_key: str):
        """Get merged categories."""
        cats_merge = set()
        decode = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)
        for storage in self.storages:
            cats = self.get_categories(storage, label_key)
            if cats is not None:
                cats = decode(cats) if isinstance(cats[0], bytes) else cats
                cats_merge.update(cats)
            else:
                codes = self.get_codes(storage, label_key)
                codes = decode(codes) if isinstance(codes[0], bytes) else codes
                cats_merge.update(codes)
        return cats_merge

    def get_categories(self, storage: StorageType, label_key: str):  # type: ignore
        """Get categories."""
        obs = storage["obs"]  # type: ignore
        if isinstance(obs, ArrayTypes):  # type: ignore
            cat_key_uns = f"{label_key}_categories"
            if cat_key_uns in storage["uns"]:  # type: ignore
                return storage["uns"][cat_key_uns]  # type: ignore
            else:
                return None
        else:
            if "__categories" in obs:
                cats = obs["__categories"]
                if label_key in cats:
                    return cats[label_key]
                else:
                    return None
            labels = obs[label_key]
            if isinstance(labels, GroupTypes):  # type: ignore
                if "categories" in labels:
                    return labels["categories"]
                else:
                    return None
            else:
                if "categories" in labels.attrs:
                    return labels.attrs["categories"]
                else:
                    return None
        return None

    def get_codes(self, storage: StorageType, label_key: str):  # type: ignore
        """Get codes."""
        obs = storage["obs"]  # type: ignore
        if isinstance(obs, ArrayTypes):  # type: ignore
            label = obs[label_key]
        else:
            label = obs[label_key]
            if isinstance(label, ArrayTypes):  # type: ignore
                return label[...]
            else:
                return label["codes"][...]

    def close(self):
        """Close connection to array streaming backend."""
        for storage in self.storages:
            if hasattr(storage, "close"):
                storage.close()
        for conn in self.conns:
            if hasattr(conn, "close"):
                conn.close()
