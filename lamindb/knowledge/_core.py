import bionty as bt


class Species(bt.Species):
    """Bionty Species.

    See `bionty.Species <https://lamin.ai/docs/bionty/bionty.Species>`__.
    """

    def __init__(self, id="name") -> None:
        super().__init__(id=id)


class Gene(bt.Gene):
    """Bionty Gene.

    See `bionty.Gene <https://lamin.ai/docs/bionty/bionty.Gene>`__.
    """

    def __init__(self, species="human", id="ensembl_gene_id") -> None:
        super().__init__(species=species, id=id)


class Protein(bt.Protein):
    """Bionty Protein.

    See `bionty.Protein <https://lamin.ai/docs/bionty/bionty.Protein>`__.
    """

    def __init__(self, species="human", id="uniprotkb_id") -> None:
        super().__init__(species=species, id=id)


class CellMarker(bt.CellMarker):
    """Bionty CellMarker.

    See `bionty.CellMarker <https://lamin.ai/docs/bionty/bionty.CellMarker>`__.
    """

    def __init__(self, species="human", id="name") -> None:
        super().__init__(species=species, id=id)


class CellType(bt.CellType):
    """Bionty CellType.

    See `bionty.CellType <https://lamin.ai/docs/bionty/bionty.CellType>`__.
    """

    def __init__(self, id="ontology_id") -> None:
        super().__init__(id=id)


class Tissue(bt.Tissue):
    """Bionty Tissue.

    See `bionty.Tissue <https://lamin.ai/docs/bionty/bionty.Tissue>`__.
    """

    def __init__(self, id="ontology_id") -> None:
        super().__init__(id=id)


class Disease(bt.Disease):
    """Bionty Disease.

    See `bionty.Disease <https://lamin.ai/docs/bionty/bionty.Disease>`__.
    """

    def __init__(self, id="ontology_id") -> None:
        super().__init__(id=id)
