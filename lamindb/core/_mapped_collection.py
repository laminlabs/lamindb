from __future__ import annotations

from collections import Counter
from functools import reduce
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Literal, Optional, Union

import numpy as np
import pandas as pd
from lamin_utils import logger
from lamindb_setup.core.upath import UPath

from .storage._backed_access import (
    ArrayType,
    ArrayTypes,
    GroupType,
    GroupTypes,
    StorageType,
    _safer_read_index,
    registry,
)

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr


class _Connect:
    def __init__(self, storage):
        if isinstance(storage, UPath):
            self.conn, self.store = registry.open("h5py", storage)
            self.to_close = True
        else:
            self.conn, self.store = None, storage
            self.to_close = False

    def __enter__(self):
        return self.store

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if not self.to_close:
            return
        if hasattr(self.store, "close"):
            self.store.close()
        if hasattr(self.conn, "close"):
            self.conn.close()


class MappedCollection:
    """Map-style collection for use in data loaders.

    This class virtually concatenates `AnnData` arrays as a `pytorch map-style dataset
    <https://pytorch.org/docs/stable/data.html#map-style-datasets>`__.

    If your `AnnData` collection is in the cloud, move them into a local cache
    first for faster access.

    `__getitem__` of the `MappedCollection` object takes a single integer index
    and returns a dictionary with the observation data sample for this index from
    the `AnnData` objects in `path_list`. The dictionary has keys for `layers_keys`
    (`.X` is in `"X"`), `obs_keys`, `obsm_keys` (under `f"obsm_{key}"`) and also `"_store_idx"`
    for the index of the `AnnData` object containing this observation sample.

    .. note::

        For a guide, see :doc:`docs:scrna5`.

        For more convenient use within :class:`~lamindb.core.MappedCollection`,
        see :meth:`~lamindb.Collection.mapped`.

        This currently only works for collections of `AnnData` objects.

        The implementation was influenced by the `SCimilarity
        <https://github.com/Genentech/scimilarity>`__ data loader.


    Args:
        path_list: A list of paths to `AnnData` objects stored in `.h5ad` or `.zarr` formats.
        layers_keys: Keys from the ``.layers`` slot. ``layers_keys=None`` or ``"X"`` in the list
            retrieves ``.X``.
        obsm_keys: Keys from the ``.obsm`` slots.
        obs_keys: Keys from the ``.obs`` slots.
        join: `"inner"` or `"outer"` virtual joins. If ``None`` is passed,
            does not join.
        encode_labels: Encode labels into integers.
            Can be a list with elements from ``obs_keys``.
        unknown_label: Encode this label to -1.
            Can be a dictionary with keys from ``obs_keys`` if ``encode_labels=True``
            or from ``encode_labels`` if it is a list.
        cache_categories: Enable caching categories of ``obs_keys`` for faster access.
        parallel: Enable sampling with multiple processes.
        dtype: Convert numpy arrays from ``.X``, ``.layers`` and ``.obsm``
    """

    def __init__(
        self,
        path_list: list[UPathStr],
        layers_keys: str | list[str] | None = None,
        obs_keys: str | list[str] | None = None,
        obsm_keys: str | list[str] | None = None,
        join: Literal["inner", "outer"] | None = "inner",
        encode_labels: bool | list[str] = True,
        unknown_label: str | dict[str, str] | None = None,
        cache_categories: bool = True,
        parallel: bool = False,
        dtype: str | None = None,
    ):
        assert join in {None, "inner", "outer"}

        if layers_keys is None:
            self.layers_keys = ["X"]
        else:
            self.layers_keys = (
                [layers_keys] if isinstance(layers_keys, str) else layers_keys
            )

        obsm_keys = [obsm_keys] if isinstance(obsm_keys, str) else obsm_keys
        self.obsm_keys = obsm_keys

        obs_keys = [obs_keys] if isinstance(obs_keys, str) else obs_keys
        self.obs_keys = obs_keys

        if isinstance(encode_labels, list):
            if len(encode_labels) == 0:
                encode_labels = False
            elif obs_keys is None or not all(
                enc_label in obs_keys for enc_label in encode_labels
            ):
                raise ValueError(
                    "All elements of `encode_labels` should be in `obs_keys`."
                )
        else:
            if encode_labels:
                encode_labels = obs_keys if obs_keys is not None else False
        self.encode_labels = encode_labels

        if encode_labels and isinstance(unknown_label, dict):
            if not all(unkey in encode_labels for unkey in unknown_label):  # type: ignore
                raise ValueError(
                    "All keys of `unknown_label` should be in `encode_labels` and `obs_keys`."
                )
        self.unknown_label = unknown_label

        self.storages = []  # type: ignore
        self.conns = []  # type: ignore
        self.parallel = parallel
        self._path_list = path_list
        self._make_connections(path_list, parallel)

        self.n_obs_list = []
        for storage in self.storages:
            with _Connect(storage) as store:
                X = store["X"]
                if isinstance(X, ArrayTypes):  # type: ignore
                    self.n_obs_list.append(X.shape[0])
                else:
                    self.n_obs_list.append(X.attrs["shape"][0])
        self.n_obs = sum(self.n_obs_list)

        self.indices = np.hstack([np.arange(n_obs) for n_obs in self.n_obs_list])
        self.storage_idx = np.repeat(np.arange(len(self.storages)), self.n_obs_list)

        self.join_vars = join
        self.var_indices = None
        self.var_joint = None
        self.n_vars_list = None
        self.n_vars = None
        if self.join_vars is not None:
            self._make_join_vars()
            self.n_vars = len(self.var_joint)

        if self.obs_keys is not None:
            if cache_categories:
                self._cache_categories(self.obs_keys)
            else:
                self._cache_cats: dict = {}
            self.encoders: dict = {}
            if self.encode_labels:
                self._make_encoders(self.encode_labels)  # type: ignore

        self._dtype = dtype
        self._closed = False

    def _make_connections(self, path_list: list, parallel: bool):
        for path in path_list:
            path = UPath(path)
            if path.exists() and path.is_file():  # type: ignore
                if parallel:
                    conn, storage = None, path
                else:
                    conn, storage = registry.open("h5py", path)
            else:
                conn, storage = registry.open("zarr", path)
            self.conns.append(conn)
            self.storages.append(storage)

    def _cache_categories(self, obs_keys: list):
        self._cache_cats = {}
        decode = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)
        for label in obs_keys:
            self._cache_cats[label] = []
            for storage in self.storages:
                with _Connect(storage) as store:
                    cats = self._get_categories(store, label)
                    if cats is not None:
                        cats = decode(cats) if isinstance(cats[0], bytes) else cats[...]
                    self._cache_cats[label].append(cats)

    def _make_encoders(self, encode_labels: list):
        for label in encode_labels:
            cats = self.get_merged_categories(label)
            encoder = {}
            if isinstance(self.unknown_label, dict):
                unknown_label = self.unknown_label.get(label, None)
            else:
                unknown_label = self.unknown_label
            if unknown_label is not None and unknown_label in cats:
                cats.remove(unknown_label)
                encoder[unknown_label] = -1
            encoder.update({cat: i for i, cat in enumerate(cats)})
            self.encoders[label] = encoder

    def _make_join_vars(self):
        var_list = []
        self.n_vars_list = []
        for storage in self.storages:
            with _Connect(storage) as store:
                vars = _safer_read_index(store["var"])
                var_list.append(vars)
                self.n_vars_list.append(len(vars))

        vars_eq = all(var_list[0].equals(vrs) for vrs in var_list[1:])
        if vars_eq:
            self.join_vars = None
            self.var_joint = var_list[0]
            return

        if self.join_vars == "inner":
            self.var_joint = reduce(pd.Index.intersection, var_list)
            if len(self.var_joint) == 0:
                raise ValueError(
                    "The provided AnnData objects don't have shared varibales.\n"
                    "Use join='outer'."
                )
            self.var_indices = [vrs.get_indexer(self.var_joint) for vrs in var_list]
        elif self.join_vars == "outer":
            self.var_joint = reduce(pd.Index.union, var_list)
            self.var_indices = [self.var_joint.get_indexer(vrs) for vrs in var_list]

    def __len__(self):
        return self.n_obs

    @property
    def shape(self):
        """Shape of the (virtually aligned) dataset."""
        return (self.n_obs, self.n_vars)

    @property
    def original_shapes(self):
        """Shapes of the underlying AnnData objects."""
        if self.n_vars_list is None:
            n_vars_list = [None] * len(self.n_obs_list)
        else:
            n_vars_list = self.n_vars_list
        return list(zip(self.n_obs_list, n_vars_list))

    def __getitem__(self, idx: int):
        obs_idx = self.indices[idx]
        storage_idx = self.storage_idx[idx]
        if self.var_indices is not None:
            var_idxs_join = self.var_indices[storage_idx]
        else:
            var_idxs_join = None

        with _Connect(self.storages[storage_idx]) as store:
            out = {}
            for layers_key in self.layers_keys:
                lazy_data = (
                    store["X"] if layers_key == "X" else store["layers"][layers_key]
                )
                out[layers_key] = self._get_data_idx(
                    lazy_data, obs_idx, self.join_vars, var_idxs_join, self.n_vars
                )
            if self.obsm_keys is not None:
                for obsm_key in self.obsm_keys:
                    lazy_data = store["obsm"][obsm_key]
                    out[f"obsm_{obsm_key}"] = self._get_data_idx(lazy_data, obs_idx)
            out["_store_idx"] = storage_idx
            if self.obs_keys is not None:
                for label in self.obs_keys:
                    if label in self._cache_cats:
                        cats = self._cache_cats[label][storage_idx]
                        if cats is None:
                            cats = []
                    else:
                        cats = None
                    label_idx = self._get_obs_idx(store, obs_idx, label, cats)
                    if label in self.encoders:
                        label_idx = self.encoders[label][label_idx]
                    out[label] = label_idx
        return out

    def _get_data_idx(
        self,
        lazy_data: ArrayType | GroupType,  # type: ignore
        idx: int,
        join_vars: Literal["inner", "outer"] | None = None,
        var_idxs_join: list | None = None,
        n_vars_out: int | None = None,
    ):
        """Get the index for the data."""
        if isinstance(lazy_data, ArrayTypes):  # type: ignore
            lazy_data_idx = lazy_data[idx]  # type: ignore
            if join_vars is None:
                result = lazy_data_idx
                if self._dtype is not None:
                    result = result.astype(self._dtype, copy=False)
            elif join_vars == "outer":
                dtype = lazy_data_idx.dtype if self._dtype is None else self._dtype
                result = np.zeros(n_vars_out, dtype=dtype)
                result[var_idxs_join] = lazy_data_idx
            else:  # inner join
                result = lazy_data_idx[var_idxs_join]
                if self._dtype is not None:
                    result = result.astype(self._dtype, copy=False)
            return result
        else:  # assume csr_matrix here
            data = lazy_data["data"]  # type: ignore
            indices = lazy_data["indices"]  # type: ignore
            indptr = lazy_data["indptr"]  # type: ignore
            s = slice(*(indptr[idx : idx + 2]))
            data_s = data[s]
            dtype = data_s.dtype if self._dtype is None else self._dtype
            if join_vars == "outer":
                lazy_data_idx = np.zeros(n_vars_out, dtype=dtype)
                lazy_data_idx[var_idxs_join[indices[s]]] = data_s
            else:
                lazy_data_idx = np.zeros(lazy_data.attrs["shape"][1], dtype=dtype)  # type: ignore
                lazy_data_idx[indices[s]] = data_s
                if join_vars == "inner":
                    lazy_data_idx = lazy_data_idx[var_idxs_join]
            return lazy_data_idx

    def _get_obs_idx(
        self,
        storage: StorageType,
        idx: int,
        label_key: str,
        categories: list | None = None,
    ):
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
        if categories is not None:
            cats = categories
        else:
            cats = self._get_categories(storage, label_key)
        if cats is not None and len(cats) > 0:
            label = cats[label]
        if isinstance(label, bytes):
            label = label.decode("utf-8")
        return label

    def get_label_weights(self, obs_keys: str | list[str]):
        """Get all weights for the given label keys."""
        if isinstance(obs_keys, str):
            obs_keys = [obs_keys]
        labels_list = []
        for label_key in obs_keys:
            labels_to_str = self.get_merged_labels(label_key).astype(str).astype("O")
            labels_list.append(labels_to_str)
        if len(labels_list) > 1:
            labels = reduce(lambda a, b: a + b, labels_list)
        else:
            labels = labels_list[0]
        labels = self.get_merged_labels(label_key)
        counter = Counter(labels)  # type: ignore
        weights = 1.0 / np.array([counter[label] for label in labels])
        return weights

    def get_merged_labels(self, label_key: str):
        """Get merged labels for `label_key` from all `.obs`."""
        labels_merge = []
        decode = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)
        for i, storage in enumerate(self.storages):
            with _Connect(storage) as store:
                codes = self._get_codes(store, label_key)
                labels = decode(codes) if isinstance(codes[0], bytes) else codes
                if label_key in self._cache_cats:
                    cats = self._cache_cats[label_key][i]
                else:
                    cats = self._get_categories(store, label_key)
                if cats is not None:
                    cats = decode(cats) if isinstance(cats[0], bytes) else cats
                    labels = cats[labels]
                labels_merge.append(labels)
        return np.hstack(labels_merge)

    def get_merged_categories(self, label_key: str):
        """Get merged categories for `label_key` from all `.obs`."""
        cats_merge = set()
        decode = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)
        for i, storage in enumerate(self.storages):
            with _Connect(storage) as store:
                if label_key in self._cache_cats:
                    cats = self._cache_cats[label_key][i]
                else:
                    cats = self._get_categories(store, label_key)
                if cats is not None:
                    cats = decode(cats) if isinstance(cats[0], bytes) else cats
                    cats_merge.update(cats)
                else:
                    codes = self._get_codes(store, label_key)
                    codes = decode(codes) if isinstance(codes[0], bytes) else codes
                    cats_merge.update(codes)
        return cats_merge

    def _get_categories(self, storage: StorageType, label_key: str):  # type: ignore
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

    def _get_codes(self, storage: StorageType, label_key: str):  # type: ignore
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
        """Close connections to array streaming backend.

        No effect if `parallel=True`.
        """
        for storage in self.storages:
            if hasattr(storage, "close"):
                storage.close()
        for conn in self.conns:
            if hasattr(conn, "close"):
                conn.close()
        self._closed = True

    @property
    def closed(self):
        """Check if connections to array streaming backend are closed.

        Does not matter if `parallel=True`.
        """
        return self._closed

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @staticmethod
    def torch_worker_init_fn(worker_id):
        """`worker_init_fn` for `torch.utils.data.DataLoader`.

        Improves performance for `num_workers > 1`.
        """
        from torch.utils.data import get_worker_info

        mapped = get_worker_info().dataset
        mapped.parallel = False
        mapped.storages = []
        mapped.conns = []
        mapped._make_connections(mapped._path_list, parallel=False)
