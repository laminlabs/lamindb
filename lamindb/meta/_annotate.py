from typing import Optional  # noqa

from bioreader import lookup
from tabulate import tabulate  # type: ignore

from .._logger import colors, logger
from ..dev.db import insert
from ..do._query import query
from ..do._update import update


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
            featureset_id=geneset_id,
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
        dobject_id,
        readout_type: lookup.READOUT_TYPES,
        readout_platform: Optional[lookup.READOUT_PLATFORMS],
    ):
        # register the readout if not yet in the database
        readout_results = query.readout_type(
            name=readout_type, platform=readout_platform
        )
        if len(readout_results) == 0:
            readout_type_id = insert.readout_type(
                name=readout_type, platform=readout_platform
            )
            logger.success(f"Registered readout_type: {readout_type_id}")
        else:
            readout_type_id = readout_results[0].id

        # query biometa associated with a dobject
        dobject_biometa = query.dobject_biometa(dobject_id=dobject_id)
        if len(dobject_biometa) > 0:
            biometa_ids = [i.biometa_id for i in dobject_biometa]
        else:
            biometa_ids = [insert.biometa(dobject_id=dobject_id)]
            logger.warning(
                f"No biometa found for dobject {dobject_id}, created biometa"
                f" {biometa_ids[0]}"
            )

        # fill in biometa entries with readout_type_id
        for biometa_id in biometa_ids:
            update.biometa(biometa_id, readout_type_id=readout_type_id)

        logger.success(
            f"{colors.blue(f'readout_type_id {readout_type_id}')} has been added to"
            f" {colors.purple(f'biometa entries {biometa_ids}')} associated with"
            f" {colors.green(f'dobject {dobject_id}')}."
        )
