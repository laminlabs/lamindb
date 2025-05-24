from __future__ import annotations

import inspect
from functools import cached_property
from itertools import chain
from typing import TYPE_CHECKING, Callable, Literal, Union

import h5py
import numpy as np
import pandas as pd
from anndata import AnnData
from anndata import __version__ as anndata_version
from anndata._core.index import _normalize_indices
from anndata._core.views import _resolve_idx
from anndata._io.h5ad import read_dataframe_legacy as read_dataframe_legacy_h5
from anndata._io.specs.registry import get_spec, read_elem, read_elem_partial
from anndata.compat import _read_attr
from fsspec.implementations.local import LocalFileSystem
from fsspec.utils import infer_compression
from lamin_utils import logger
from lamindb_setup.core.upath import infer_filesystem
from packaging import version
from upath import UPath

if TYPE_CHECKING:
    from collections.abc import Mapping

    from fsspec.core import OpenFile
    from lamindb_setup.core.types import UPathStr


anndata_version_parse = version.parse(anndata_version)

if anndata_version_parse < version.parse("0.9.0"):
    from anndata._core.index import Index
else:
    from anndata.compat import Index

if anndata_version_parse < version.parse("0.10.0"):
    if anndata_version_parse < version.parse("0.9.1"):
        logger.warning(
            "Full backed capabilities are not available for this version of anndata,"
            " please install anndata>=0.9.1."
        )

    from anndata._core.sparse_dataset import SparseDataset

    # try csr for groups with no encoding_type
    class CSRDataset(SparseDataset):
        @property
        def format_str(self) -> str:
            return "csr"

    def sparse_dataset(group):
        return SparseDataset(group)

else:
    if anndata_version_parse >= version.parse("0.11.0"):
        from anndata._core.sparse_dataset import (  # type: ignore
            _CSRDataset as CSRDataset,
        )
    else:
        from anndata._core.sparse_dataset import CSRDataset  # type: ignore
    from anndata._core.sparse_dataset import (
        BaseCompressedSparseDataset as SparseDataset,
    )
    from anndata._core.sparse_dataset import sparse_dataset  # type: ignore

    def _check_group_format(*args):
        pass

    CSRDataset._check_group_format = _check_group_format


# zarr and CSRDataset have problems with full selection
def _subset_sparse(sparse_ds: CSRDataset | SparseDataset, indices):
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


class AccessRegistry:
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
registry = AccessRegistry()


@registry.register_open("h5py")
def open(filepath: UPathStr, mode: str = "r", compression: str | None = "infer"):
    fs, file_path_str = infer_filesystem(filepath)
    # we don't open compressed files directly because we need fsspec to uncompress on .open
    compression = (
        infer_compression(file_path_str) if compression == "infer" else compression
    )
    if isinstance(fs, LocalFileSystem) and compression is None:
        assert mode in {"r", "r+", "a", "w", "w-"}, f"Unknown mode {mode}!"  #  noqa: S101
        return None, h5py.File(file_path_str, mode=mode)
    if mode == "r":
        conn_mode = "rb"
    elif mode == "w":
        conn_mode = "wb"
    elif mode == "a":
        conn_mode = "ab"
    else:
        raise ValueError(f"Unknown mode {mode}! Should be 'r', 'w' or 'a'.")
    conn = fs.open(file_path_str, mode=conn_mode, compression=compression)
    try:
        storage = h5py.File(conn, mode=mode)
    except Exception as e:
        conn.close()
        raise e
    return conn, storage


@registry.register("h5py")
def read_dataframe(elem: h5py.Dataset | h5py.Group):
    if isinstance(elem, h5py.Dataset):
        return read_dataframe_legacy_h5(elem)
    else:
        return read_elem(elem)


@registry.register("h5py")
def safer_read_partial(elem, indices):
    is_dataset = isinstance(elem, h5py.Dataset)
    indices_inverse: list | None = None
    encoding_type = get_spec(elem).encoding_type
    # h5py selection for datasets requires sorted indices
    if is_dataset or encoding_type == "dataframe":
        indices_increasing = []
        indices_inverse = []
        for indices_dim in indices:
            # should be integer or bool
            # ignore bool or increasing unique integers
            if (
                isinstance(indices_dim, np.ndarray)
                and indices_dim.dtype != "bool"
                and not np.all(np.diff(indices_dim) > 0)
            ):
                idx_unique, idx_inverse = np.unique(indices_dim, return_inverse=True)
                indices_increasing.append(idx_unique)
                indices_inverse.append(idx_inverse)
            else:
                indices_increasing.append(indices_dim)
                indices_inverse.append(None)
        indices = tuple(indices_increasing)
        if all(idx is None for idx in indices_inverse):
            indices_inverse = None
    result = None
    if encoding_type == "":
        if is_dataset:
            dims = len(elem.shape)
            if dims == 2:
                result = elem[indices]
            elif dims == 1:
                if indices[0] == slice(None):
                    result = elem[indices[1]]
                elif indices[1] == slice(None):
                    result = elem[indices[0]]
        elif isinstance(elem, h5py.Group):
            try:
                ds = CSRDataset(elem)
                result = _subset_sparse(ds, indices)
            except Exception as e:
                logger.debug(
                    f"Encountered an exception while attempting to subset a sparse dataset by indices.\n{e}"
                )
        if result is None:
            raise ValueError(
                "Can not get a subset of the element of type"
                f" {type(elem).__name__} with an empty spec."
            )
    else:
        result = read_elem_partial(elem, indices=indices)
    if indices_inverse is None:
        return result
    else:
        if indices_inverse[0] is None:
            if len(result.shape) == 2:
                return result[:, indices_inverse[1]]
            else:
                return result[indices_inverse[1]]
        elif indices_inverse[1] is None:
            if isinstance(result, pd.DataFrame):
                return result.iloc[indices_inverse[0]]
            else:
                return result[indices_inverse[0]]
        else:
            return result[tuple(indices_inverse)]


@registry.register("h5py")
def keys(storage: h5py.File):
    attrs_keys: dict[str, list] = {}
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

    from ._zarr import get_zarr_store

    ArrayTypes.append(zarr.Array)
    GroupTypes.append(zarr.Group)
    StorageTypes.append(zarr.Group)

    @registry.register_open("zarr")
    def open(filepath: UPathStr, mode: Literal["r", "r+", "a", "w", "w-"] = "r"):
        assert mode in {"r", "r+", "a", "w", "w-"}, f"Unknown mode {mode}!"  #  noqa: S101

        store = get_zarr_store(filepath)
        storage = zarr.open(store, mode=mode)
        conn = None
        return conn, storage

    @registry.register("zarr")
    def read_dataframe(elem: Union[zarr.Array, zarr.Group]):  # noqa
        if isinstance(elem, zarr.Array):
            return read_dataframe_legacy_zarr(elem)
        else:
            return read_elem(elem)

    @registry.register("zarr")
    def safer_read_partial(elem, indices):
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
                except Exception as e:
                    logger.debug(
                        f"Encountered an exception while attempting to subset a sparse dataset by indices.\n{e}"
                    )
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
    def keys(storage: zarr.Group):
        if hasattr(storage, "_sync_iter"):  # zarr v3
            paths = storage._sync_iter(storage.store.list())
        else:
            paths = storage.store.keys()  # zarr v2

        attrs_keys: dict[str, list] = {}
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
            obj = registry.safer_read_partial(self.storage["obs"], indices=indices)  # type: ignore
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
            obj = registry.safer_read_partial(self.storage["var"], indices=indices)  # type: ignore
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
            if isinstance(oidx, np.ndarray) or oidx != slice(None):
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
            if isinstance(var_raw, ArrayTypes):
                attrs_keys["var"] = list(var_raw.dtype.fields.keys())
            else:
                # for some reason list(var_raw.keys()) is very slow for zarr
                # maybe also directly get keys from the underlying mapper
                attrs_keys["var"] = list(var_raw)
            if "varm" in storage_raw:
                varm_keys_raw = list(storage_raw["varm"])
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
        connection: OpenFile | None,
        storage: StorageType,
        filename: str,
    ):
        self._conn = connection
        self.storage = storage

        self._attrs_keys = registry.keys(self.storage)

        self._name = filename

        self._obs_names = _safer_read_index(self.storage["obs"])  # type: ignore
        self._var_names = _safer_read_index(self.storage["var"])  # type: ignore

        self._closed = False

    def close(self):
        """Closes the connection."""
        if hasattr(self, "storage") and hasattr(self.storage, "close"):
            self.storage.close()
        if hasattr(self, "_conn") and hasattr(self._conn, "close"):
            self._conn.close()
        self._closed = True

    @property
    def closed(self):
        return self._closed

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
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


# get the number of observations in an anndata object or file fast and safely
def _anndata_n_observations(object: UPathStr | AnnData) -> int | None:
    if isinstance(object, AnnData):
        return object.n_obs

    try:
        objectpath = UPath(object)
        suffix = objectpath.suffix
        conn_module = {".h5ad": "h5py", ".zarr": "zarr"}.get(suffix, suffix[1:])
        conn, storage = registry.open(conn_module, objectpath, mode="r")
    except Exception as e:
        logger.warning(f"Could not open {object} to read n_observations: {e}")
        return None

    n_observations: int | None = None
    try:
        obs = storage["obs"]
        if isinstance(obs, GroupTypes):  # type: ignore
            if "_index" in obs.attrs:
                elem_key = _read_attr(obs.attrs, "_index")
            else:
                elem_key = next(iter(obs))
            elem = obs[elem_key]
            if isinstance(elem, ArrayTypes):  # type: ignore
                n_observations = elem.shape[0]
            else:
                # assume standard obs group
                n_observations = elem["codes"].shape[0]
        else:
            n_observations = obs.shape[0]
    except Exception as e:
        logger.warning(f"Could not read n_observations from anndata {object}: {e}")
    finally:
        if hasattr(storage, "close"):
            storage.close()
        if hasattr(conn, "close"):
            conn.close()
    return n_observations
