import pandas as pd
import pytest
from lnschema_bionty import CellType

import lamindb as ln


@pytest.fixture(scope="module")
def instance():
    ln.setup.init(storage="test_parse", schema="bionty")
    yield
    ln.setup.delete("test_parse")


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


def test_parse_name(instance, df):
    result = ln.parse(df.cell_type, CellType.name)
    ids = [i.ontology_id for i in result]
    assert len(result) == 3
    assert set(ids) == set(["CL:0000182", None, "CL:0000084"])


def test_parse_ontology_id(instance, df):
    result = ln.parse(df.cell_type_id, CellType.ontology_id)
    names = [i.name for i in result]
    assert len(result) == 2
    assert set(names) == set(["T cell", "hepatocyte"])


def test_parse_df(instance, df):
    result = ln.parse(
        df, {"cell_type": CellType.name, "cell_type_id": CellType.ontology_id}
    )
    names = [i.name for i in result]
    ids = [i.ontology_id for i in result]
    extras = [i.definition for i in result]
    assert len(result) == 3
    assert set(names) == set(["T cell", "hepatocyte", "my new cell type"])
    # converts empty string to None
    assert set(ids) == set(["CL:0000084", "CL:0000182", None])
    # additional fields are not populated
    assert set(extras) == set([None, None, None])
