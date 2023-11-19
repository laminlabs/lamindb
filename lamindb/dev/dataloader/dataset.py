from collections import Counter
from os import PathLike
from pathlib import Path
from typing import List, Optional, Union

import numpy as np

from ..storage._backed_access import ArrayTypes, GroupTypes, StorageType, registry


class ListDataset:
    def __init__(
        self,
        pth_list: List[Union[str, PathLike]],
        labels: Optional[Union[str, List[str]]] = None,
        encode_labels: bool = True,
    ):
        self.storages = []
        self.conns = []
        for pth in pth_list:
            pth = Path(pth) if not isinstance(pth, Path) else pth
            if pth.exists() and pth.is_file():
                conn, storage = registry.open("h5py", pth)
            else:
                conn, storage = registry.open("zarr", pth)
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
        if isinstance(labels, str):
            labels = [labels]
        self.labels = labels
        if self.labels is not None:
            if self.encode_labels:
                self.encoders = []
                for label in self.labels:
                    cats = self.get_merged_categories(label)
                    self.encoders.append({cat: i for i, cat in enumerate(cats)})

    def __len__(self):
        return self.n_obs

    def __getitem__(self, idx):
        obs_idx = self.indices[idx]
        storage = self.storages[self.storage_idx[idx]]
        out = [self.get_data_idx(storage, obs_idx)]
        if self.labels is not None:
            for i, label in enumerate(self.labels):
                label_idx = self.get_label_idx(storage, obs_idx, label)
                if self.encode_labels:
                    label_idx = self.encoders[i][label_idx]
                out.append(label_idx)
        return out

    def get_data_idx(
        self, storage: StorageType, idx: int, layer_key: Optional[str] = None  # type: ignore # noqa
    ):
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

    def get_labels_weights(self, label_key: str):
        labels = self.get_merged_labels(label_key)
        counter = Counter(labels)  # type: ignore
        weights = np.array([counter[label] for label in labels]) / len(labels)
        return weights

    def get_merged_labels(self, label_key: str):
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
        obs = storage["obs"]  # type: ignore
        if isinstance(obs, ArrayTypes):  # type: ignore
            label = obs[label_key]
        else:
            label = obs[label_key]
            if isinstance(label, ArrayTypes):  # type: ignore
                return label[...]
            else:
                return label["codes"][...]
