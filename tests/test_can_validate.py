import bionty as bt
import lamindb as ln
import pytest


# some validate tests are in test_queryset
def test_inspect():
    ln.FeatureSet.filter().all().delete()
    bt.Gene.filter().all().delete()
    bt.settings.organism = "human"
    result = bt.Gene.inspect("TCF7", "symbol")
    assert result.validated == []

    bt.Gene.from_source(symbol="TCF7").save()
    result = bt.Gene.inspect("TCF7")
    assert bt.Gene.validate("TCF7", organism="human")
    result = bt.Gene.inspect(["TCF7", "ABC1"], "symbol")
    assert result.validated == ["TCF7"]

    # clean up
    bt.Gene.filter().all().delete()


def test_standardize():
    bt.settings.organism = "human"

    # synonym not in the database
    result = bt.Gene.standardize(["ABC1", "PDCD1"])
    assert result == ["HEATR6", "PDCD1"]

    result = bt.Gene.standardize(["ABC1", "PDCD1"], field=bt.Gene.symbol)
    assert result == ["HEATR6", "PDCD1"]

    mapper = bt.Gene.standardize(["ABC1", "PDCD1"], return_mapper=True)
    assert mapper == {"ABC1": "HEATR6"}

    # synonym already in the database
    bt.Gene.from_source(symbol="LMNA").save()
    mapper = bt.Gene.standardize(["ABC1", "LMN1"], return_mapper=True)
    assert mapper == {"LMN1": "LMNA", "ABC1": "HEATR6"}
    assert bt.Gene.standardize(["LMNA"]) == ["LMNA"]
    assert bt.Gene.standardize("LMNA") == "LMNA"
    assert bt.Gene.standardize(["LMN1"], return_mapper=True) == {"LMN1": "LMNA"}


def test_standardize_public_aware():
    result = bt.Gene.standardize(["ABC1", "PDCD1"], public_aware=False)
    assert result == ["ABC1", "PDCD1"]


def test_add_remove_synonym():
    bt.CellType.filter().all().delete()
    # a registry that cannot validate
    source = bt.Source.filter(organism="human").first()
    with pytest.raises(AttributeError):
        source.add_synonym("syn")

    # a registry that doesn't have a synonyms column
    user = ln.User.filter(handle="testuser1").one()
    with pytest.raises(NotImplementedError):
        user.add_synonym("syn")

    cell_types = bt.CellType.from_values(["T cell", "B cell"], "name")
    ln.save(cell_types)
    tcell = bt.CellType.filter(name="T cell").one()
    bcell = bt.CellType.filter(name="B cell").one()
    tcell.add_synonym(["my cell type"])
    tcell.add_synonym("")
    tcell.add_synonym([])
    assert "my cell type" in tcell.synonyms
    with pytest.raises(ValueError):
        bcell.add_synonym("my cell type")
    with pytest.raises(ValueError):
        tcell.add_synonym("my|celltype")

    tcell.remove_synonym("my cell type")
    assert "my cell type" not in tcell.synonyms

    bcell.synonyms = None
    bcell.save()
    tcell.synonyms = None
    tcell.save()
    tcell.add_synonym("")
    tcell.add_synonym([""])
    tcell.add_synonym([])
    tcell.add_synonym(["my cell type"])
    tcell.add_synonym("")
    tcell.add_synonym([""])
    tcell.add_synonym([])
    assert tcell.synonyms == "my cell type"
    tcell.remove_synonym("my cell type")

    # clean up
    bt.CellType.filter().all().delete()


def test_set_abbr():
    bt.CellType.filter().all().delete()
    bt.CellType(name="my cell type").save()
    record = bt.CellType.filter(name="my cell type").one()
    # if abbr is name, do not add to synonyms
    record.set_abbr("my cell type")
    assert record.abbr == "my cell type"
    assert record.synonyms is None

    record.set_abbr("myct")
    assert record.abbr == "myct"
    assert "myct" in record.synonyms

    source = bt.Source.filter(organism="human").first()
    with pytest.raises(AttributeError) as error:
        source.set_abbr("abbr")
    assert (
        error.exconly() == "AttributeError: 'Source' object has no attribute 'set_abbr'"
    )

    record.delete()
