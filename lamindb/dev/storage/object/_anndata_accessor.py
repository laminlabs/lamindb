from functools import cached_property
from typing import Dict, Mapping, Union

import h5py
import pandas as pd
import zarr
from anndata._core.index import Index, _normalize_indices
from anndata._core.sparse_dataset import SparseDataset
from anndata._core.views import _resolve_idx
from anndata._io.specs.methods import read_indices
from anndata._io.specs.registry import get_spec, read_elem, read_elem_partial
from anndata.compat import _read_attr
from lamindb_setup.dev.upath import infer_filesystem as _infer_filesystem
from lnschema_core import File
from lnschema_core.models import filepath_from_file_or_folder

from ._subset_anndata import _read_dataframe


def _try_backed_full(elem):
    # think what to do for compatibility with old var and obs
    if isinstance(elem, (h5py.Dataset, zarr.Array)):
        return elem

    if isinstance(elem, (h5py.Group, zarr.Group)):
        encoding_type = get_spec(elem).encoding_type
        if encoding_type in ("csr_matrix", "csc_matrix"):
            return SparseDataset(elem)
        if "h5sparse_format" in elem.attrs:
            return SparseDataset(elem)

    return read_elem(elem)


class _MapAccessor:
    def __init__(self, elem, name, indices=None):
        self.elem = elem
        self.indices = indices
        self.name = name

    def __getitem__(self, key):
        if self.indices is None:
            return _try_backed_full(self.elem[key])
        else:
            return read_elem_partial(self.elem[key], indices=self.indices)

    def keys(self):
        return list(self.elem.keys())

    def __repr__(self):
        """Description of the _MapAccessor object."""
        descr = f"Accessor for the AnnData attribute {self.name}"
        descr += f"\n  with keys: {self.keys()}"
        return descr


class _AnnDataAttrsMixin:
    storage: Union[h5py.File, zarr.Group]
    _attrs_keys: Mapping[str, list]

    @cached_property
    def obs(self) -> pd.DataFrame:
        if "obs" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[0], slice(None))
            return read_elem_partial(self.storage["obs"], indices=indices)
        else:
            return _read_dataframe(self.storage["obs"])

    @cached_property
    def var(self) -> pd.DataFrame:
        if "var" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[1], slice(None))
            return read_elem_partial(self.storage["var"], indices=indices)
        else:
            return _read_dataframe(self.storage["var"])

    @cached_property
    def uns(self):
        if "uns" not in self._attrs_keys:
            return None
        return read_elem(self.storage["uns"])

    @cached_property
    def X(self):
        indices = getattr(self, "indices", None)
        if indices is not None:
            return read_elem_partial(self.storage["X"], indices=indices)
        else:
            return _try_backed_full(self.storage["X"])

    @cached_property
    def obsm(self):
        if "obsm" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[0], slice(None))
        return _MapAccessor(self.storage["obsm"], "obsm", indices)

    @cached_property
    def varm(self):
        if "varm" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[1], slice(None))
        return _MapAccessor(self.storage["varm"], "varm", indices)

    @cached_property
    def obsp(self):
        if "obsp" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[0], indices[0])
        return _MapAccessor(self.storage["obsp"], "obsp", indices)

    @cached_property
    def varp(self):
        if "varp" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[1], indices[1])
        return _MapAccessor(self.storage["varp"], "varp", indices)

    @cached_property
    def layers(self):
        if "layers" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        return _MapAccessor(self.storage["layers"], "layers", indices)

    @property
    def obs_names(self):
        return self._obs_names

    @property
    def var_names(self):
        return self._var_names

    @cached_property
    def shape(self):
        return len(self._obs_names), len(self._var_names)


class AnnDataAccessorSubset(_AnnDataAttrsMixin):
    def __init__(self, storage, indices, attrs_keys, obs_names, var_names, ref_shape):
        self.storage = storage
        self.indices = indices

        self._attrs_keys = attrs_keys
        self._obs_names, self._var_names = obs_names, var_names

        self._ref_shape = ref_shape

    def __getitem__(self, index: Index):
        """Access a subset of the underlying AnnData object."""
        oidx, vidx = _normalize_indices(index, self._obs_names, self._var_names)
        new_obs_names, new_var_names = self._obs_names[oidx], self._var_names[vidx]
        if self.indices is not None:
            oidx = _resolve_idx(self.indices[0], oidx, self._ref_shape[0])
            vidx = _resolve_idx(self.indices[1], vidx, self._ref_shape[1])
        return type(self)(
            self.storage,
            (oidx, vidx),
            self._attrs_keys,
            new_obs_names,
            new_var_names,
            self._ref_shape,
        )

    def __repr__(self):
        """Description of the object."""
        n_obs, n_vars = self.shape
        descr = f"{type(self).__name__} object with n_obs × n_vars = {n_obs} × {n_vars}"
        for attr, keys in self._attrs_keys.items():
            descr += f"\n  {attr}: {keys}"
        return descr

    @cached_property
    def raw(self):
        if "raw" not in self._attrs_keys:
            return None
        prepare_indices = None
        if self.indices is not None:
            oidx = self.indices[0]
            if oidx != slice(None):
                prepare_indices = oidx, slice(None)
        return AnnDataRawAccessor(
            self.storage["raw"],
            prepare_indices,
            None,
            self._obs_names,
            None,
            self._ref_shape[0],
        )


class AnnDataRawAccessor(AnnDataAccessorSubset):
    def __init__(
        self, storage_raw, indices, attrs_keys, obs_names, var_names, ref_shape
    ):
        var_raw = storage_raw["var"]

        if var_names is None:
            var_names = read_elem(var_raw[_read_attr(var_raw.attrs, "_index")])

        if isinstance(ref_shape, int):
            ref_shape = ref_shape, len(var_names)
        elif isinstance(ref_shape, tuple) and len(ref_shape) < 2:
            ref_shape = ref_shape[0], len(var_names)

        if attrs_keys is None:
            attrs_keys = {}
            if isinstance(var_raw, (h5py.Dataset, zarr.Array)):
                attrs_keys["var"] = list(var_raw.dtype.fields.keys())
            else:
                # for some reason list(var_raw.keys()) is very slow for zarr
                # maybe also directly get keys from the underlying mapper
                attrs_keys["var"] = [key for key in var_raw]
            if "varm" in storage_raw:
                varm_keys_raw = [key for key in storage_raw["varm"]]
                if len(varm_keys_raw) > 0:
                    attrs_keys["varm"] = varm_keys_raw

        super().__init__(
            storage_raw, indices, attrs_keys, obs_names, var_names, ref_shape
        )

    @property
    def raw(self):
        raise AttributeError


class AnnDataAccessor(_AnnDataAttrsMixin):
    def __init__(self, file: File):
        fs, file_path_str = _infer_filesystem(filepath_from_file_or_folder(file))

        if file.suffix == ".h5ad":
            self._conn = fs.open(file_path_str, mode="rb")
            try:
                self.storage = h5py.File(self._conn, mode="r")
            except Exception as e:
                self._conn.close()
                raise e
            self._attrs_keys = _keys_h5(self.storage)
        elif file.suffix in (".zarr", ".zrad"):
            self._conn = None
            mapper = fs.get_mapper(file_path_str, check=True)
            self.storage = zarr.open(mapper, mode="r")
            self._attrs_keys = _keys_zarr(self.storage)
        else:
            raise ValueError(
                f"file should have .h5ad, .zarr or .zrad suffix, not {file.suffix}."
            )

        self._name = file.name

        self._obs_names, self._var_names = read_indices(self.storage)

    def __del__(self):
        """Closes the connection."""
        if self._conn is not None:
            self.storage.close()
            self._conn.close()

    def __getitem__(self, index: Index) -> AnnDataAccessorSubset:
        """Access a subset of the underlying AnnData object."""
        oidx, vidx = _normalize_indices(index, self._obs_names, self._var_names)
        new_obs_names, new_var_names = self._obs_names[oidx], self._var_names[vidx]
        return AnnDataAccessorSubset(
            self.storage,
            (oidx, vidx),
            self._attrs_keys,
            new_obs_names,
            new_var_names,
            self.shape,
        )

    def __repr__(self):
        """Description of the AnnDataAccessor object."""
        n_obs, n_vars = self.shape
        descr = f"AnnDataAccessor object with n_obs × n_vars = {n_obs} × {n_vars}"
        descr += f"\n  constructed for the AnnData object {self._name}"
        for attr, keys in self._attrs_keys.items():
            descr += f"\n    {attr}: {keys}"
        return descr

    @cached_property
    def raw(self):
        if "raw" not in self._attrs_keys:
            return None
        return AnnDataRawAccessor(
            self.storage["raw"], None, None, self._obs_names, None, self.shape[0]
        )


# this is needed because accessing zarr.Group.keys() directly is very slow
def _keys_zarr(storage: zarr.Group):
    paths = storage._store.keys()

    attrs_keys: Dict[str, list] = {}
    obs_var_arrays = []

    for path in paths:
        if path in (".zattrs", ".zgroup"):
            continue
        parts = path.split("/")
        if len(parts) < 2:
            continue
        attr = parts[0]
        key = parts[1]

        if attr == "X":
            continue

        if attr in ("obs", "var"):
            if attr in obs_var_arrays:
                continue
            if key == ".zarray":
                attrs_keys.pop(attr, None)
                obs_var_arrays.append(attr)

        if attr not in attrs_keys:
            attrs_keys[attr] = []

        if key in (".zattrs", ".zgroup", ".zarray"):
            continue
        attr_keys = attrs_keys[attr]
        if key not in attr_keys:
            attr_keys.append(key)

    for attr in obs_var_arrays:
        attrs_keys[attr] = list(storage[attr].dtype.fields.keys())

    return {attr: keys for attr, keys in attrs_keys.items() if len(keys) > 0}


def _keys_h5(storage: h5py.File):
    attrs_keys: Dict[str, list] = {}
    for attr in storage.keys():
        if attr == "X":
            continue
        attr_obj = storage[attr]
        if attr in ("obs", "var") and isinstance(attr_obj, h5py.Dataset):
            keys = list(attr_obj.dtype.fields.keys())
        else:
            keys = list(attr_obj.keys())
        if len(keys) > 0:
            attrs_keys[attr] = keys
    return attrs_keys
