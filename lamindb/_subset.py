from typing import List, Optional, Union

from anndata import AnnData, concat
from lamin_logger import logger
from lnschema_core import DObject

from .dev.object import _subset_anndata_dobject
from .dev.object._lazy_field import LazySelector

SUFFIXES = (".h5ad", ".zarr")


def subset(
    dobjects: Union[List[DObject], DObject],
    query_obs: Optional[Union[List[str], str, LazySelector]] = None,
    query_var: Optional[Union[List[str], str, LazySelector]] = None,
    use_concat: bool = False,
    concat_args: Optional[dict] = None,
) -> Union[List[AnnData], AnnData, None]:
    """Subset AnnData dobjects and stream results into memory.

    Args:
        dobjects: A `DObject` or a list of `DObjects` containing `AnnData` objects
        to subset and load into memory.
        query_obs: The pandas query string to evaluate on `.obs` of each
        underlying `AnnData` object.
        query_var: The pandas query string to evaluate on `.var` of each
        underlying `AnnData` object.
        use_concat: If `True`, applies `anndata.concat` on
        the returned `AnnData` objects.
        concat_args: Arguments for concatenation.
    """
    if isinstance(dobjects, DObject):
        dobjects = [dobjects]

    n_dobjects = len(dobjects)
    if isinstance(query_obs, list):
        if len(query_obs) != n_dobjects:
            raise ValueError("query_obs list should be the same length as dobjects.")
    else:
        query_obs = [query_obs] * n_dobjects  # type: ignore
    if isinstance(query_var, list):
        if len(query_var) != n_dobjects:
            raise ValueError("query_var list should be the same length as dobjects.")
    else:
        query_var = [query_var] * n_dobjects  # type: ignore

    adatas = []
    # todo: implement parallel processing
    for i, dobject in enumerate(dobjects):
        if dobject.suffix not in SUFFIXES:
            logger.warning(f"DObject {dobject.id} is not an AnnData object, ignoring.")
            continue
        result = _subset_anndata_dobject(dobject, query_obs[i], query_var[i])
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
