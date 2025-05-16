import bionty as bt
import lamindb as ln
import pytest
from lamindb.errors import ValidationError


# some validate tests are in test_queryset
def test_inspect():
    ln.Schema.filter().delete()
    bt.Gene.filter().delete()
    bt.settings.organism = "human"
    result = bt.Gene.inspect("TCF7", "symbol")
    assert result.validated == []

    bt.Gene.from_source(symbol="TCF7").save()
    result = bt.Gene.inspect("TCF7")
    assert bt.Gene.validate("TCF7", organism="human")
    result = bt.Gene.inspect(["TCF7", "ABC1"], "symbol")
    assert result.validated == ["TCF7"]

    # clean up
    bt.Gene.filter().delete()


# if a record was added to the DB via a different source
# it will still be validated because it's in the DB
def test_inspect_source():
    source1 = bt.Source.get(entity="bionty.CellType", name="cl")
    source2 = bt.CellType.add_source(source="cl", version="2022-08-16")
    bt.CellType.from_source(name="T cell", source=source1).save()
    assert bt.CellType.inspect("T-cell", source=source2, mute=True).synonyms_mapper == {
        "T-cell": "T cell"
    }
    assert (
        bt.CellType.inspect(
            "T-cell", source=source2, mute=True, strict_source=True
        ).synonyms_mapper
        == {}
    )
    assert bt.CellType.validate("T cell", source=source2, mute=True).sum() == 1
    assert (
        bt.CellType.validate(
            "T cell", source=source2, mute=True, strict_source=True
        ).sum()
        == 0
    )
    assert bt.CellType.standardize("T-cell", source=source2, mute=True) == "T cell"
    # here still standardized because of bionty
    assert (
        bt.CellType.standardize("T-cell", source=source2, mute=True, strict_source=True)
        == "T cell"
    )
    bt.CellType.filter().delete()


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


def test_standardize_source_aware():
    result = bt.Gene.standardize(["ABC1", "PDCD1"], source_aware=False)
    assert result == ["ABC1", "PDCD1"]


def test_add_remove_synonym():
    bt.CellType.filter().delete()

    # a registry that doesn't have a synonyms column
    user = ln.User.get(handle="testuser1")
    with pytest.raises(NotImplementedError):
        user.add_synonym("syn")

    cell_types = bt.CellType.from_values(["T cell", "B cell"], "name")
    ln.save(cell_types)
    tcell = bt.CellType.get(name="T cell")
    bcell = bt.CellType.get(name="B cell")
    tcell.add_synonym(["my cell type"])
    tcell.add_synonym("")
    tcell.add_synonym([])
    assert "my cell type" in tcell.synonyms
    with pytest.raises(ValidationError):
        bcell.add_synonym("my cell type")
    with pytest.raises(ValidationError):
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
    bt.CellType.filter().delete()


def test_set_abbr():
    bt.CellType.filter().delete()
    bt.CellType(name="my cell type").save()
    record = bt.CellType.get(name="my cell type")
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


def test_validate_int():
    result = ln.User.validate([1, 2], field=ln.User.id)
    assert result.sum() == 1


def test_synonym_mapping():
    # only name field can be standardized
    bt.Gene.from_source(symbol="TNFRSF4", organism="human").save()

    result = bt.Gene.inspect(
        ["CD134", "TNFRSF4"], field=bt.Gene.symbol, organism="human"
    )
    assert result.synonyms_mapper == {"CD134": "TNFRSF4"}

    result = bt.Gene.inspect(
        ["CD134", "TNFRSF4"], field=bt.Gene.ensembl_gene_id, organism="human"
    )
    assert result.synonyms_mapper == {}

    bt.Gene.filter().delete()
