import bionty as bt


def Gene(species="human", id="ensembl_gene_id"):
    """Returns `bionty.Gene <https://lamin.ai/docs/bionty/bionty.Gene>`__."""
    return bt.Gene(species=species, id=id)


def Protein(species="human", id="uniprotkb_id"):
    """Returns `bionty.Protein <https://lamin.ai/docs/bionty/bionty.Protein>`__."""
    return bt.Protein(species=species, id=id)


# fmt: off
def CellMarker(species="human", id="name"):
    """Returns `bionty.CellMarker <https://lamin.ai/docs/bionty/bionty.CellMarker>`__."""  # noqa
    return bt.CellMarker(species=species, id=id)
# fmt: on
