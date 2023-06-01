from typing import List, Optional, Union

from anndata import AnnData, concat
from lamin_logger import logger
from lnschema_core import File

from .object import LazySelector, _subset_anndata_file

SUFFIXES = (".h5ad", ".zarr")


def subset(
    files: Union[List[File], File],
    query_obs: Optional[Union[List[str], str, LazySelector]] = None,
    query_var: Optional[Union[List[str], str, LazySelector]] = None,
    use_concat: bool = False,
    concat_args: Optional[dict] = None,
) -> Union[List[AnnData], AnnData, None]:
    """Subset AnnData files and stream results into memory.

    Args:
        files: A `File` or a list of `Files` containing `AnnData` objects
        to subset and load into memory.
        query_obs: The pandas query string to evaluate on `.obs` of each
        underlying `AnnData` object.
        query_var: The pandas query string to evaluate on `.var` of each
        underlying `AnnData` object.
        use_concat: If `True`, applies `anndata.concat` on
        the returned `AnnData` objects.
        concat_args: Arguments for concatenation.
    """
    if isinstance(files, File):
        files = [files]

    n_files = len(files)
    if isinstance(query_obs, list):
        if len(query_obs) != n_files:
            raise ValueError("query_obs list should be the same length as files.")
    else:
        query_obs = [query_obs] * n_files  # type: ignore
    if isinstance(query_var, list):
        if len(query_var) != n_files:
            raise ValueError("query_var list should be the same length as files.")
    else:
        query_var = [query_var] * n_files  # type: ignore

    adatas = []
    # todo: implement parallel processing
    for i, file in enumerate(files):
        if file.suffix not in SUFFIXES:
            logger.warning(f"File {file.id} is not an AnnData object, ignoring.")
            continue
        result = _subset_anndata_file(file, query_obs[i], query_var[i])
        if result is not None:
            adatas.append(result)

    if not use_concat:
        return adatas
    # concat branch here
    n_adatas = len(adatas)
    if n_adatas == 0:
        return None
    elif n_adatas == 1:
        return adatas[0]
    else:
        concat_args = {} if concat_args is None else concat_args
        return concat(adatas, **concat_args)
