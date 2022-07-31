from typing import Optional  # noqa

from bioreader import vocabulary as vc
from tabulate import tabulate  # type: ignore

from .._logger import colors, logger
from ..dev.db import insert
from ..do._query import query


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
    def gene(
        cls,
        dobject_id,
        values: dict,
        species: str,
        geneset_name: str = None,
    ):
        """Annotate genes."""
        geneset_id = insert.genes(
            genes_dict=values, geneset_name=geneset_name, species=species
        )

        # use the geneset_id and readout_type_id to create an entry in biometa
        biometa_id = insert.biometa(
            dobject_id=dobject_id,
            geneset_id=geneset_id,
        )

        logs = [[str(geneset_id), str(biometa_id)]]
        log_table = tabulate(
            logs,
            headers=[
                colors.green("geneset.id"),
                colors.purple("biometa.id"),
            ],
            tablefmt="pretty",
        )
        logger.success(
            f"Annotated data {dobject_id} with the following features:\n{log_table}",
        )

    @classmethod
    def readout_type(
        cls,
        readout_type: vc.READOUT_TYPES,
        readout_platform: Optional[vc.READOUT_PLATFORMS],
    ):
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

        # query biometa

        # fill in biometa entries with readout_type_id

        return readout_type_id
