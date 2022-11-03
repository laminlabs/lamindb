import bioreadout
from lamin_logger import colors, logger

from lamindb.dev.db._add import add
from lamindb.dev.db._select import select
from lamindb.schema import wetlab


class link:
    """Link metadata.

    Guide: :doc:`/db/guide/link`.
    """

    @classmethod
    def readout(cls, dobject_id, efo_id: str):
        """Link readout."""
        readout = select(wetlab.readout, efo_id=efo_id).one_or_none()
        if readout is None:
            assert sum(i.isdigit() for i in efo_id) == 7
            readout = add(wetlab.readout(**bioreadout.readout(efo_id=efo_id)))

        # select biometa associated with a dobject
        dobject_biometas = select(wetlab.dobject_biometa, dobject_id=dobject_id).all()
        if len(dobject_biometas) > 0:
            biometa_ids = [i.biometa_id for i in dobject_biometas]
            biometas = (
                select(wetlab.biometa).where(wetlab.biometa.id.in_(biometa_ids)).all()
            )
            for biometa in biometas:
                biometa.readout_id = readout.id
            add(biometas)
        else:
            add(wetlab.biometa(readout_id=readout.id))
            logger.warning(
                f"No biometa found for dobject {dobject_id}, created biometa"
                f" {biometa_ids[0]}"
            )

        logger.success(
            f"Added {colors.blue(f'readout_id {readout.id}')} to"
            f" {colors.purple(f'biometa {biometa_ids}')} linked to"
            f" {colors.green(f'dobject {dobject_id}')}."
        )

    @classmethod
    def biometa(cls, dobject_id: str, biometa_id: str):
        """Link a dobject to a biometa."""
        dobject_biometas = select(
            wetlab.dobject_biometa, dobject_id=dobject_id, biometa_id=biometa_id
        ).all()
        if len(dobject_biometas) > 0:
            raise AssertionError(
                f"dobject {dobject_id} is already linked to biometa {biometa_id}!"
            )
        else:
            add(wetlab.dobject_biometa(dobject_id=dobject_id, biometa_id=biometa_id))
