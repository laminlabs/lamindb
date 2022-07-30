from typing import Literal, Optional  # noqa

from bioreader import vocabulary as vc
from tabulate import tabulate  # type: ignore

from .._logger import colors, logger
from ..dev.db import insert
from ..dev.file import h5ad_to_anndata
from ..do import query
from ..schema import core


def anndata_to_df(adata, obs_or_var):
    """Get a df from AnnData."""
    if obs_or_var == "var":
        df = adata.var
    elif obs_or_var == "obs":
        df = adata.obs
    else:
        raise KeyError("Please specify obs or var!")
    return df


class annotate:
    """Feature annotation."""

    @classmethod
    def readouts(
        cls,
        dobject: core.dobject,
        species: str,
        readout_entity: str,
        readout_type: vc.READOUT_TYPES,
        readout_platform: Optional[vc.READOUT_PLATFORMS],
        column=None,
        obs_or_var=None,
        readout_set_name: str = None,
    ):
        """Annotate genes."""
        filekey = f"{dobject.id}-{dobject.v}{dobject.file_suffix}"

        if dobject.file_suffix == ".h5ad":
            adata = h5ad_to_anndata(filekey)
            df = anndata_to_df(adata, obs_or_var=obs_or_var)

            # create a geneset entry
            readouts = (
                df.index.unique().values if column is None else df[column].unique()
            )
            readout_set_id = insert.readouts(
                readouts=readouts,
                readout_entity=readout_entity,
                readout_set_name=readout_set_name,
                species=species,
            )

            # register the readout if not yet in the database
            readout_results = query.readout_type(
                name=readout_type, platform=readout_platform
            )
            if len(readout_results) == 0:
                readout_type_id = insert.readout_type(
                    name=readout_type, platform=readout_platform
                )
            else:
                readout_type_id = readout_results[0].id

            # use the geneset_id and readout_type_id to create an entry in biometa
            biometa_id = insert.biometa(
                dobject_id=dobject.id,
                readout_type_id=readout_type_id,
                geneset_id=readout_set_id,
            )

            logs = [[str(readout_set_id), str(readout_type_id), str(biometa_id)]]
            log_table = tabulate(
                logs,
                headers=[
                    colors.green("geneset.id"),
                    colors.blue("readout_type.id"),
                    colors.purple("biometa.id"),
                ],
                tablefmt="pretty",
            )
            logger.success(
                f"{colors.bold('Annotated the following features')}:\n{log_table}",
            )
        else:
            raise NotImplementedError
