import bionty as bt
import lamindb as ln
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame(
        (
            ["T cell", "CL:0000084"],
            ["hepatocyte", "CL:0000182"],
            ["my new cell type", ""],
        ),
        columns=["cell_type", "cell_type_id"],
    )


def test_from_values_name(df):
    bt.CellType.filter().delete()
    assert df["cell_type"].tolist() == ["T cell", "hepatocyte", "my new cell type"]
    # create records from bionty
    result = bt.CellType.from_values(df.cell_type, "name")
    ids = [i.ontology_id for i in result]
    assert len(result) == 2
    assert set(ids) == {"CL:0000084", "CL:0000182"}
    assert result[0].public_source.entity == "CellType"

    # wrong field type
    with pytest.raises(TypeError):
        result = bt.CellType.from_values(df.cell_type, field=bt.CellType)


def test_from_values_ontology_id(df):
    assert df["cell_type_id"].tolist() == ["CL:0000084", "CL:0000182", ""]
    result = bt.CellType.from_values(df.cell_type_id, "ontology_id")
    names = {i.name for i in result}
    assert len(result) == 2
    assert names == {"T cell", "hepatocyte"}
    assert result[0].public_source.entity == "CellType"


def test_from_values_multiple_match():
    records = bt.Gene.from_values(["ABC1", "PDCD1"], bt.Gene.symbol, organism="human")
    assert len(records) == 3


def test_from_values_organism():
    from lnschema_bionty import Gene, settings

    settings._organism = None

    with pytest.raises(AssertionError):
        Gene.from_values(["ABC1"], Gene.symbol)

    settings.organism = "human"
    values = ["ABC1"]
    curated_values = Gene.public().standardize(values)
    records = Gene.from_values(curated_values, Gene.symbol)
    assert records[0].ensembl_gene_id == "ENSG00000068097"

    settings.organism = "mouse"
    values = ["ABC1"]
    curated_values = Gene.public().standardize(values)
    records = Gene.from_values(curated_values, Gene.symbol)
    assert records[0].ensembl_gene_id == "ENSMUSG00000015243"


def test_get_or_create_records():
    names = ["ulabel" + str(i) for i in range(25)]
    labels = [ln.ULabel(name=name) for name in names]
    ln.save(labels)
    # more than 20 existing values
    labels = ln.ULabel.from_values(names, field="name")
    assert len(labels) == 25


def test_from_values_synonyms_aware():
    bt.CellType.from_public(name="T cell").save(parents=False)
    # existing validated values
    records = bt.CellType.from_values(["T cell"], "name")
    assert len(records) == 1
    assert records[0].name == "T cell"
    assert isinstance(records[0].public_source, bt.PublicSource)
    # existing validated values and synonyms
    records = bt.CellType.from_values(["T cell", "T-cell"], "name")
    assert len(records) == 1
    assert records[0].name == "T cell"
    assert isinstance(records[0].public_source, bt.PublicSource)
    # bionty values and synonyms
    records = bt.CellType.from_values(["B-cell", "B cell"], "name")
    assert len(records) == 1
    assert records[0].name == "B cell"
    assert isinstance(records[0].public_source, bt.PublicSource)
    # all possibilities of validated values
    records = bt.CellType.from_values(
        ["T cell", "T-cell", "t cell", "B cell", "B-cell"], "name"
    )
    assert len(records) == 2
    names = [r.name for r in records]
    assert set(names) == {"T cell", "B cell"}
    assert isinstance(records[0].public_source, bt.PublicSource)
    assert isinstance(records[1].public_source, bt.PublicSource)
    # non-validated values
    records = bt.CellType.from_values(["T cell", "mycell"], "name")
    assert len(records) == 1
    assert records[0].name == "T cell"
    assert isinstance(records[0].public_source, bt.PublicSource)
    assert records[0].ontology_id == "CL:0000084"
    bt.CellType.filter().all().delete()
