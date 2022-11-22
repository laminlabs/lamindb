from typing import List, Optional, Union

from anndata import AnnData, concat
from lnschema_core import DObject

from .dev.object import _subset_dobject


def subset(
    dobjects: List[DObject],
    query_obs: Optional[str] = None,
    query_var: Optional[str] = None,
    use_concat: bool = False,
    concat_args: Optional[dict] = None,
) -> Union[List[AnnData], AnnData, None]:
    adatas = []
    # todo: implement parallel processing
    for dobject in dobjects:
        result = _subset_dobject(dobject, query_obs, query_var)
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
