import lnschema_bionty as lb
import pandas as pd
import pytest

import lamindb as ln  # noqa


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
    assert df["cell_type_id"].tolist() == ["CL:0000084", "CL:0000182", ""]
    result = lb.CellType.from_values(df.cell_type, "name")
    ids = [i.ontology_id for i in result]
    assert len(result) == 3
    assert set(ids) == {"CL:0000084", "CL:0000182", None}
    for r in result:
        if r.ontology_id == "CL:0000084":
            assert r.bionty_source.entity == "CellType"
        if r.ontology_id is None:
            assert r.bionty_source is None
    with pytest.raises(TypeError):
        result = lb.CellType.from_values(df.cell_type, field=lb.CellType)


def test_from_values_ontology_id(df):
    result = lb.CellType.from_values(df.cell_type_id, "ontology_id")
    names = {i.name for i in result}
    assert len(result) == 2
    assert names == {"T cell", "hepatocyte"}
    assert result[0].bionty_source.entity == "CellType"


def test_from_values_multiple_match():
    records = lb.Gene.from_values(["ABC1", "PDCD1"], lb.Gene.symbol, species="human")
    assert len(records) == 3


def test_from_values_species():
    from lnschema_bionty import Gene, settings

    settings._species = None

    with pytest.raises(AssertionError):
        Gene.from_values(["ABC1"], Gene.symbol)

    settings.species = "human"
    values = ["ABC1"]
    curated_values = Gene.bionty().map_synonyms(values)
    print(curated_values)
    records = Gene.from_values(curated_values, Gene.symbol)
    assert records[0].ensembl_gene_id == "ENSG00000068097"

    settings.species = "mouse"
    values = ["ABC1"]
    curated_values = Gene.bionty().map_synonyms(values)
    print(curated_values)
    records = Gene.from_values(curated_values, Gene.symbol)
    assert records[0].ensembl_gene_id == "ENSMUSG00000015243"


def test_get_or_create_records():
    labels = ln.Label.from_values(["label" + str(i) for i in range(25)], field="name")
    ln.save(labels)
    # more than 20 existing values
    ln.Label.from_values(["label" + str(i) for i in range(25)], field="name")
