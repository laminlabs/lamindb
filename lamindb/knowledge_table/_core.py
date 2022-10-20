import bionty as bt


def Gene(species="human", id="ensembl_gene_id"):
    """Knowledge table of genes."""
    return bt.Gene(species=species, id=id)


def Protein(species="human", id="uniprotkb_id"):
    """Knowledge table of proteins."""
    return bt.Protein(species=species, id=id)


def CellMarker(species="human", id="name"):
    """Knowledge table of cell markers."""
    return bt.CellMarker(species=species, id=id)
