from typing import Optional, Union

import h5py
import zarr
from anndata import AnnData
from anndata._io.h5ad import read_dataframe_legacy as read_dataframe_legacy_h5
from anndata._io.specs.methods import _read_partial
from anndata._io.specs.registry import read_elem, read_elem_partial
from anndata._io.zarr import read_dataframe_legacy as read_dataframe_legacy_zarr
from lamindb_setup.dev.upath import infer_filesystem as _infer_filesystem
from lnschema_core.models import File, filepath_from_file_or_folder

from ._lazy_field import LazySelector


def _read_dataframe(elem: Union[zarr.Array, h5py.Dataset, zarr.Group, h5py.Group]):
    if isinstance(elem, zarr.Array):
        return read_dataframe_legacy_zarr(elem)
    elif isinstance(elem, h5py.Dataset):
        return read_dataframe_legacy_h5(elem)
    else:
        return read_elem(elem)


def _indices(base_indices, select_indices):
    if len(base_indices) == len(select_indices):
        return slice(None)
    else:
        return list(base_indices.get_indexer(select_indices))


def _subset_adata_storage(
    storage: Union[zarr.Group, h5py.File],
    query_obs: Optional[Union[str, LazySelector]] = None,
    query_var: Optional[Union[str, LazySelector]] = None,
) -> Union[AnnData, None]:
    with storage as access:
        obs = _read_dataframe(access["obs"])
        var = _read_dataframe(access["var"])

        if query_obs is not None:
            if hasattr(query_obs, "evaluate"):
                obs_result = obs[query_obs.evaluate(obj=obs)]  # type: ignore
            else:
                obs_result = obs.query(query_obs)
        else:
            obs_result = obs

        if query_var is not None:
            if hasattr(query_var, "evaluate"):
                var_result = var[query_var.evaluate(obj=var)]  # type: ignore
            else:
                var_result = var.query(query_var)
        else:
            var_result = var

        if obs_result.index.empty or var_result.index.empty:
            return None

        obs_idx = _indices(obs.index, obs_result.index)
        var_idx = _indices(var.index, var_result.index)

        prepare_adata = {}
        prepare_adata["obs"] = obs_result
        prepare_adata["var"] = var_result
        X = read_elem_partial(access["X"], indices=(obs_idx, var_idx))
        prepare_adata["X"] = X
        if "obsm" in access:
            obsm = _read_partial(access["obsm"], indices=(obs_idx, slice(None)))
            prepare_adata["obsm"] = obsm
        if "varm" in access:
            varm = _read_partial(access["varm"], indices=(var_idx, slice(None)))
            prepare_adata["varm"] = varm
        if "obsp" in access:
            obsp = _read_partial(access["obsp"], indices=(obs_idx, obs_idx))
            prepare_adata["obsp"] = obsp
        if "varp" in access:
            varp = _read_partial(access["varp"], indices=(var_idx, var_idx))
            prepare_adata["varp"] = varp
        if "layers" in access:
            layers = _read_partial(access["layers"], indices=(obs_idx, var_idx))
            prepare_adata["layers"] = layers
        if "uns" in access:
            prepare_adata["uns"] = read_elem(access["uns"])

        return AnnData(**prepare_adata)


def _subset_anndata_file(
    file: File,
    query_obs: Optional[Union[str, LazySelector]] = None,
    query_var: Optional[Union[str, LazySelector]] = None,
) -> Union[AnnData, None]:
    file_path = filepath_from_file_or_folder(file)
    fs, file_path_str = _infer_filesystem(file_path)

    adata: Union[AnnData, None] = None

    if file.suffix == ".h5ad":
        with fs.open(file_path_str, mode="rb") as file:
            storage = h5py.File(file, mode="r")
            adata = _subset_adata_storage(storage, query_obs, query_var)
    elif file.suffix == ".zarr":
        mapper = fs.get_mapper(file_path_str, check=True)
        storage = zarr.open(mapper, mode="r")
        adata = _subset_adata_storage(storage, query_obs, query_var)

    return adata
