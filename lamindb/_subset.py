from typing import List, Optional, Union

from anndata import AnnData, concat
from lamin_logger import logger
from lnschema_core import DObject

from .dev.object import _subset_anndata_dobject

SUFFIXES = (".h5ad", ".zarr")


def subset(
    dobjects: Union[List[DObject], DObject],
    query_obs: Optional[str] = None,
    query_var: Optional[str] = None,
    use_concat: bool = False,
    concat_args: Optional[dict] = None,
) -> Union[List[AnnData], AnnData, None]:
    """Subset AnnData dobjects and stream results into memory."""
    if isinstance(dobjects, DObject):
        dobjects = [dobjects]

    adatas = []
    # todo: implement parallel processing
    for dobject in dobjects:
        if dobject.suffix not in SUFFIXES:
            logger.warning(f"DObject {dobject.id} is not an AnnData object, ignoring.")
            continue
        result = _subset_anndata_dobject(dobject, query_obs, query_var)
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
