import bionty as bt


# fmt: off
def Gene(species="human", id="ensembl_gene_id"):
    """Returns a `bionty.Gene instance <https://lamin.ai/docs/bionty/bionty.Gene>`__."""  # noqa
    return bt.Gene(species=species, id=id)


def Protein(species="human", id="uniprotkb_id"):
    """Returns a `bionty.Protein instance <https://lamin.ai/docs/bionty/bionty.Protein>`__."""  # noqa
    return bt.Protein(species=species, id=id)


def CellMarker(species="human", id="name"):
    """Returns a `bionty.CellMarker instance <https://lamin.ai/docs/bionty/bionty.CellMarker>`__."""  # noqa
    return bt.CellMarker(species=species, id=id)
# fmt: on
