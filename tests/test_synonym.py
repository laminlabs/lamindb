import lnschema_bionty as lb
import pytest

import lamindb as ln  # noqa


def test_map_synonyms():
    lb.settings.species = "human"

    # synonym not in the database
    result = lb.Gene.map_synonyms(["ABC1", "PDCD1"])
    assert result == ["HEATR6", "PDCD1"]

    result = lb.Gene.map_synonyms(["ABC1", "PDCD1"], field=lb.Gene.symbol)
    assert result == ["HEATR6", "PDCD1"]

    mapper = lb.Gene.map_synonyms(["ABC1", "PDCD1"], return_mapper=True)
    assert mapper == {"ABC1": "HEATR6"}

    # synonym already in the database
    lb.Gene.from_bionty(symbol="LMNA").save()
    mapper = lb.Gene.map_synonyms(["ABC1", "LMN1"], return_mapper=True)
    assert mapper == {"LMN1": "LMNA", "ABC1": "HEATR6"}
    assert lb.Gene.map_synonyms(["LMNA"]) == ["LMNA"]


def test_add_remove_synonym():
    bionty_source = lb.BiontySource.filter(species="human").first()
    with pytest.raises(NotImplementedError):
        bionty_source.add_synonym("syn")

    cell_types = lb.CellType.from_values(["T cell", "B cell"], "name")
    ln.save(cell_types)
    tcell = lb.CellType.filter(name="T cell").one()
    bcell = lb.CellType.filter(name="B cell").one()
    tcell.add_synonym(["my cell type"])
    tcell.add_synonym("")
    tcell.add_synonym([])
    assert "my cell type" in tcell.synonyms
    with pytest.raises(SystemExit) as excinfo:
        bcell.add_synonym("my cell type")
    assert excinfo.value.code == AssertionError
    with pytest.raises(AssertionError):
        tcell.add_synonym("my|celltype")

    tcell.remove_synonym("my cell type")
    assert "my cell type" not in tcell.synonyms

    tcell.synonyms = None
    tcell.save()
    tcell.add_synonym(["my cell type"])
    tcell.synonyms == "my cell type"
