import inspect
from dataclasses import dataclass
from functools import cached_property
from itertools import chain
from pathlib import Path
from typing import Callable, Dict, Mapping, Union

import h5py
import numpy as np
import pandas as pd
from anndata import AnnData
from anndata import __version__ as anndata_version
from anndata._core.index import Index, _normalize_indices
from anndata._core.views import _resolve_idx
from anndata._io.h5ad import read_dataframe_legacy as read_dataframe_legacy_h5
from anndata._io.specs.registry import get_spec, read_elem, read_elem_partial
from anndata.compat import _read_attr
from fsspec.core import OpenFile
from lamindb_setup.dev.upath import UPath, infer_filesystem
from lnschema_core import File
from packaging import version

from lamindb.dev.storage.file import filepath_from_file

if version.parse(anndata_version) < version.parse("0.10.0"):
    from anndata._core.sparse_dataset import SparseDataset

    # try csr for groups with no encoding_type
    class CSRDataset(SparseDataset):
        @property
        def format_str(self) -> str:
            return "csr"

    def sparse_dataset(group):
        return SparseDataset(group)

else:
    from anndata._core.sparse_dataset import (
        BaseCompressedSparseDataset as SparseDataset,
    )
    from anndata._core.sparse_dataset import CSRDataset, sparse_dataset  # type: ignore

    def _check_group_format(*args):
        pass

    CSRDataset._check_group_format = _check_group_format


# zarr and CSRDataset have problems with full selection
def _subset_sparse(sparse_ds: Union[CSRDataset, SparseDataset], indices):
    has_arrays = isinstance(indices[0], np.ndarray) or isinstance(
        indices[1], np.ndarray
    )
    if not has_arrays and indices == (slice(None), slice(None)):
        return sparse_ds.to_memory()
    else:
        return sparse_ds[indices]


def get_module_name(obj):
    return inspect.getmodule(obj).__name__.partition(".")[0]


def _records_to_df(obj):
    if isinstance(obj, pd.DataFrame):
        return obj

    if hasattr(obj, "dtype") and obj.dtype.names is not None:
        formats = []
        for name, (dt, _) in obj.dtype.fields.items():
            if dt.char == "S":
                new_dt = str(dt).replace("S", "U")
            else:
                new_dt = dt
            formats.append((name, new_dt))
        df = pd.DataFrame(obj.astype(formats, copy=False))
        for index_name in ("index", "_index"):
            if index_name in df.columns:
                return df.set_index(index_name)
            return df
    else:
        return obj


class Registry:
    def __init__(self):
        self._registry = {}
        self._openers = {}

    def register_open(self, module: str):
        def wrapper(func: Callable):
            self._openers[module] = func
            return func

        return wrapper

    def open(self, module: str, *args, **kwargs):
        if module in self._openers:
            return self._openers[module](*args, **kwargs)
        else:
            raise ValueError(f"Module {module} not found, please install it.")

    def register(self, module: str):
        def wrapper(func: Callable):
            func_name = func.__name__
            if func_name not in self._registry:
                self._registry[func_name] = {}
            self._registry[func_name][module] = func
            return func

        return wrapper

    def __getattr__(self, func_name: str):
        def wrapper(*args, **kwargs):
            func_registry = self._registry[func_name]
            for arg in chain(args, kwargs.values()):
                arg_module = get_module_name(arg)
                if arg_module in func_registry:
                    return func_registry[arg_module](*args, **kwargs)
            raise ValueError(f"{func_name} is not registered for this module.")

        return wrapper


# storage specific functions should be registered and called through the registry
registry = Registry()


@registry.register_open("h5py")
def open(filepath: Union[UPath, Path, str]):
    fs, file_path_str = infer_filesystem(filepath)
    conn = fs.open(file_path_str, mode="rb")
    try:
        storage = h5py.File(conn, mode="r")
    except Exception as e:
        conn.close()
        raise e
    return conn, storage


@registry.register("h5py")
def read_dataframe(elem: Union[h5py.Dataset, h5py.Group]):
    if isinstance(elem, h5py.Dataset):
        return read_dataframe_legacy_h5(elem)
    else:
        return read_elem(elem)


@registry.register("h5py")
def safer_read_partial(elem, indices):
    if get_spec(elem).encoding_type == "":
        if isinstance(elem, h5py.Dataset):
            dims = len(elem.shape)
            if dims == 2:
                return elem[indices]
            elif dims == 1:
                if indices[0] == slice(None):
                    return elem[indices[1]]
                elif indices[1] == slice(None):
                    return elem[indices[0]]
        elif isinstance(elem, h5py.Group):
            try:
                ds = CSRDataset(elem)
                return _subset_sparse(ds, indices)
            except Exception:
                pass
        raise ValueError(
            "Can not get a subset of the element of type"
            f" {type(elem).__name__} with an empty spec."
        )
    else:
        return read_elem_partial(elem, indices=indices)


@registry.register("h5py")
def keys(storage: h5py.File):
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


ArrayTypes = [h5py.Dataset]
GroupTypes = [h5py.Group]
StorageTypes = [h5py.File]


ZARR_INSTALLED = False
try:
    import zarr

    ZARR_INSTALLED = True
except ImportError:
    pass

if ZARR_INSTALLED:
    from anndata._io.zarr import read_dataframe_legacy as read_dataframe_legacy_zarr

    ArrayTypes.append(zarr.Array)
    GroupTypes.append(zarr.Group)
    StorageTypes.append(zarr.Group)

    @registry.register_open("zarr")
    def open(filepath: Union[UPath, Path, str]):  # noqa
        fs, file_path_str = infer_filesystem(filepath)
        conn = None
        storage = zarr.open(fs.get_mapper(file_path_str, check=True), mode="r")
        return conn, storage

    @registry.register("zarr")
    def read_dataframe(elem: Union[zarr.Array, zarr.Group]):  # noqa
        if isinstance(elem, zarr.Array):
            return read_dataframe_legacy_zarr(elem)
        else:
            return read_elem(elem)

    @registry.register("zarr")
    def safer_read_partial(elem, indices):  # noqa
        encoding_type = get_spec(elem).encoding_type
        if encoding_type == "":
            if isinstance(elem, zarr.Array):
                dims = len(elem.shape)
                if dims == 2:
                    return elem.oindex[indices]
                elif dims == 1:
                    if indices[0] == slice(None):
                        return elem.oindex[indices[1]]
                    elif indices[1] == slice(None):
                        return elem.oindex[indices[0]]
            elif isinstance(elem, zarr.Group):
                try:
                    ds = CSRDataset(elem)
                    return _subset_sparse(ds, indices)
                except Exception:
                    pass
            raise ValueError(
                "Can not get a subset of the element of type"
                f" {type(elem).__name__} with an empty spec."
            )
        else:
            if encoding_type in ("csr_matrix", "csc_matrix"):
                ds = sparse_dataset(elem)
                return _subset_sparse(ds, indices)
            else:
                return read_elem_partial(elem, indices=indices)

    # this is needed because accessing zarr.Group.keys() directly is very slow
    @registry.register("zarr")
    def keys(storage: zarr.Group):  # noqa
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


ArrayTypes = tuple(ArrayTypes)  # type: ignore
GroupTypes = tuple(GroupTypes)  # type: ignore
StorageTypes = tuple(StorageTypes)  # type: ignore


ArrayType = Union[ArrayTypes]  # type: ignore
GroupType = Union[GroupTypes]  # type: ignore
StorageType = Union[StorageTypes]  # type: ignore


def _to_memory(elem):
    if isinstance(elem, ArrayTypes):
        return elem[()]
    elif isinstance(elem, SparseDataset):
        return elem.to_memory()
    else:
        return elem


def _try_backed_full(elem):
    # think what to do for compatibility with old var and obs
    if isinstance(elem, ArrayTypes):
        return elem

    if isinstance(elem, GroupTypes):
        encoding_type = get_spec(elem).encoding_type
        if encoding_type in ("csr_matrix", "csc_matrix"):
            return sparse_dataset(elem)
        if "h5sparse_format" in elem.attrs:
            return sparse_dataset(elem)
        if encoding_type == "" and "indptr" in elem:
            return CSRDataset(elem)

    return read_elem(elem)


def _safer_read_index(elem):
    if isinstance(elem, GroupTypes):
        return pd.Index(read_elem(elem[_read_attr(elem.attrs, "_index")]))
    elif isinstance(elem, ArrayTypes):
        indices = None
        for index_name in ("index", "_index"):
            if index_name in elem.dtype.names:
                indices = elem[index_name]
                break
        if indices is not None and len(indices) > 0:
            if isinstance(indices[0], bytes):
                indices = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)(indices)
            return pd.Index(indices)
        else:
            raise ValueError("Indices not found.")
    else:
        raise ValueError(f"Unknown elem type {type(elem)} when reading indices.")


class _MapAccessor:
    def __init__(self, elem, name, indices=None):
        self.elem = elem
        self.indices = indices
        self.name = name

    def __getitem__(self, key):
        if self.indices is None:
            return _try_backed_full(self.elem[key])
        else:
            return registry.safer_read_partial(self.elem[key], indices=self.indices)

    def keys(self):
        return list(self.elem.keys())

    def __repr__(self):
        """Description of the _MapAccessor object."""
        descr = f"Accessor for the AnnData attribute {self.name}"
        descr += f"\n  with keys: {self.keys()}"
        return descr


class _AnnDataAttrsMixin:
    storage: StorageType
    _attrs_keys: Mapping[str, list]

    @cached_property
    def obs(self) -> pd.DataFrame:
        if "obs" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[0], slice(None))
            obj = registry.safer_read_partial(self.storage["obs"], indices=indices)  # type: ignore # noqa
            return _records_to_df(obj)
        else:
            return registry.read_dataframe(self.storage["obs"])  # type: ignore

    @cached_property
    def var(self) -> pd.DataFrame:
        if "var" not in self._attrs_keys:
            return None
        indices = getattr(self, "indices", None)
        if indices is not None:
            indices = (indices[1], slice(None))
            obj = registry.safer_read_partial(self.storage["var"], indices=indices)  # type: ignore # noqa
            return _records_to_df(obj)
        else:
            return registry.read_dataframe(self.storage["var"])  # type: ignore

    @cached_property
    def uns(self):
        if "uns" not in self._attrs_keys:
            return None
        return read_elem(self.storage["uns"])

    @cached_property
    def X(self):
        indices = getattr(self, "indices", None)
        if indices is not None:
            return registry.safer_read_partial(self.storage["X"], indices=indices)
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

    def to_dict(self):
        prepare_adata = {}

        prepare_adata["X"] = _to_memory(self.X)

        if "uns" in self._attrs_keys:
            prepare_adata["uns"] = self.uns

        for attr in ("obs", "var"):
            if attr in self._attrs_keys:
                prepare_adata[attr] = getattr(self, attr)

        for attr in ("obsm", "varm", "obsp", "varp", "layers"):
            if attr in self._attrs_keys:
                prepare_adata[attr] = {}
                get_attr = getattr(self, attr)
                for key in self._attrs_keys[attr]:
                    prepare_adata[attr][key] = _to_memory(get_attr[key])

        if "raw" in self._attrs_keys:
            prepare_adata["raw"] = self.raw.to_dict()

        return prepare_adata

    def to_memory(self):
        adata = AnnData(**self.to_dict())
        return adata


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
            var_names = _safer_read_index(var_raw)

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
    """Cloud-backed AnnData."""

    def __init__(
        self,
        connection: Union[OpenFile, None],
        storage: StorageType,
        filename: str,
    ):
        self._conn = connection
        self.storage = storage

        self._attrs_keys = registry.keys(self.storage)

        self._name = filename

        self._obs_names = _safer_read_index(self.storage["obs"])  # type: ignore
        self._var_names = _safer_read_index(self.storage["var"])  # type: ignore

    def close(self):
        """Closes the connection."""
        if self._conn is not None:
            self.storage.close()
            self._conn.close()

    def __del__(self):
        """Closes the connection on deletion."""
        self.close()

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


@dataclass
class BackedAccessor:
    """h5py.File or zarr.Group accessor."""

    connection: OpenFile
    """The connection."""
    storage: StorageType
    """The storage access."""


def backed_access(
    file_or_filepath: Union[File, Path]
) -> Union[AnnDataAccessor, BackedAccessor]:
    if isinstance(file_or_filepath, File):
        filepath = filepath_from_file(file_or_filepath)
    else:
        filepath = file_or_filepath
    name = filepath.name

    if filepath.suffix in (".h5", ".hdf5", ".h5ad"):
        conn, storage = registry.open("h5py", filepath)
    elif filepath.suffix in (".zarr", ".zrad"):
        conn, storage = registry.open("zarr", filepath)
    else:
        raise ValueError(
            "file should have .h5, .hdf5, .h5ad, .zarr or .zrad suffix, not"
            f" {filepath.suffix}."
        )

    if filepath.suffix in (".h5ad", ".zrad"):
        return AnnDataAccessor(conn, storage, name)
    else:
        if get_spec(storage).encoding_type == "anndata":
            return AnnDataAccessor(conn, storage, name)
        else:
            return BackedAccessor(conn, storage)
