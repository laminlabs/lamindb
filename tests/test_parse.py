import pandas as pd
import pytest

import lamindb as ln


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


def test_parse_name(df):
    from lnschema_bionty import CellType

    result = ln.parse(df.cell_type, CellType.name)
    ids = [i.ontology_id for i in result]
    assert len(result) == 3
    assert ids == ["CL:0000084", "CL:0000182", None]
    assert result[0].bionty_source.entity == "CellType"
    assert result[2].bionty_source is None


def test_parse_ontology_id(df):
    from lnschema_bionty import CellType

    result = ln.parse(df.cell_type_id, CellType.ontology_id)
    names = [i.name for i in result]
    assert len(result) == 2
    assert names == ["T cell", "hepatocyte"]
    assert result[0].bionty_source.entity == "CellType"


def test_parse_df(df):
    from lnschema_bionty import CellType

    result = ln.parse(
        df, {"cell_type": CellType.name, "cell_type_id": CellType.ontology_id}
    )
    names = [i.name for i in result]
    ids = [i.ontology_id for i in result]
    extras = [i.definition for i in result]
    assert len(result) == 3
    assert names == ["T cell", "hepatocyte", "my new cell type"]
    # converts empty string to None
    assert ids == ["CL:0000084", "CL:0000182", None]
    # additional fields are not populated
    assert extras == [None, None, None]
