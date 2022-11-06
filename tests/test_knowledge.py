from lamindb.knowledge import CellMarker, Gene, Protein, Species


def test_import():
    Species(id="common_name")
    Gene()
    Protein()
    CellMarker()
