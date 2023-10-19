import lnschema_bionty as lb
import pytest

import lamindb as ln  # noqa


# some validate tests are in test_queryset
def test_inspect():
    lb.Gene.filter().all().delete()
    lb.settings.organism = "human"
    result = lb.Gene.inspect("TCF7", "symbol")
    assert result.validated == []

    lb.Gene.from_bionty(symbol="TCF7").save()
    result = lb.Gene.inspect("TCF7")
    assert lb.Gene.validate("TCF7", organism="human")
    result = lb.Gene.inspect(["TCF7", "ABC1"], "symbol")
    assert result.validated == ["TCF7"]

    # clean up
    lb.Gene.filter().all().delete()


def test_standardize():
    lb.settings.organism = "human"

    # synonym not in the database
    result = lb.Gene.standardize(["ABC1", "PDCD1"])
    assert result == ["HEATR6", "PDCD1"]

    result = lb.Gene.standardize(["ABC1", "PDCD1"], field=lb.Gene.symbol)
    assert result == ["HEATR6", "PDCD1"]

    mapper = lb.Gene.standardize(["ABC1", "PDCD1"], return_mapper=True)
    assert mapper == {"ABC1": "HEATR6"}

    # synonym already in the database
    lb.Gene.from_bionty(symbol="LMNA").save()
    mapper = lb.Gene.standardize(["ABC1", "LMN1"], return_mapper=True)
    assert mapper == {"LMN1": "LMNA", "ABC1": "HEATR6"}
    assert lb.Gene.standardize(["LMNA"]) == ["LMNA"]
    assert lb.Gene.standardize("LMNA") == "LMNA"
    assert lb.Gene.standardize(["LMN1"], return_mapper=True) == {"LMN1": "LMNA"}


def test_standardize_bionty_aware():
    result = lb.Gene.standardize(["ABC1", "PDCD1"], bionty_aware=False)
    assert result == ["ABC1", "PDCD1"]


def test_add_remove_synonym():
    lb.CellType.filter().all().delete()
    # a registry that cannot validate
    bionty_source = lb.BiontySource.filter(organism="human").first()
    with pytest.raises(AttributeError):
        bionty_source.add_synonym("syn")

    # a registry that doesn't have a synonyms column
    user = ln.User.filter(handle="testuser1").one()
    with pytest.raises(NotImplementedError):
        user.add_synonym("syn")

    cell_types = lb.CellType.from_values(["T cell", "B cell"], "name")
    ln.save(cell_types, parents=False)
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
    tcell.synonyms == "my cell type"
    tcell.remove_synonym("my cell type")

    # clean up
    lb.CellType.filter().all().delete()


def test_set_abbr():
    lb.CellType.filter().all().delete()
    lb.CellType(name="my cell type").save(parents=False)
    record = lb.CellType.filter(name="my cell type").one()
    # if abbr is name, do not add to synonyms
    record.set_abbr("my cell type")
    assert record.abbr == "my cell type"
    assert record.synonyms is None

    record.set_abbr("myct")
    assert record.abbr == "myct"
    assert "myct" in record.synonyms

    bionty_source = lb.BiontySource.filter(organism="human").first()
    with pytest.raises(AttributeError) as error:
        bionty_source.set_abbr("abbr")
    assert (
        error.exconly()
        == "AttributeError: 'BiontySource' object has no attribute 'set_abbr'"
    )  # noqa

    record.delete()
