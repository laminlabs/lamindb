from lamindb.knowledge import CellMarker, Gene, Protein, Species


def test_import():
    Species(id="name")
    Gene()
    Protein()
    CellMarker()
