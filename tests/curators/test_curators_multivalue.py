import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.core.exceptions import ValidationError


@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "sample_id": [["sample1", "sample2"], ["sample2"], ["sample3"]],
            "dose": [[1.2, 2.3], [1.2], [2.3]],
            "cell_type": [["B cell", "T cell"], ["B cell"], ["T cell"]],
            "tissue": [["blood", "pulmo"], ["blood"], ["lung"]],
        }
    )


@pytest.fixture(scope="module")
def lists_schema():
    schema = ln.Schema(
        name="lists schema cat",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="dose", dtype=float).save(),
            ln.Feature(name="cell_type", dtype=str).save(),
            ln.Feature(name="tissue", dtype=bt.Tissue).save(),
        ],
    ).save()

    yield schema

    schema.delete()
    ln.Feature.filter().delete()
    bt.Tissue.filter().delete()


def test_curator_df_multivalue(df, lists_schema):
    curator = ln.curators.DataFrameCurator(df, lists_schema)
    with pytest.raises(ValidationError):
        curator.validate()
    assert curator._cat_manager._cat_vectors["tissue"]._validated == ["blood", "lung"]
    assert curator._cat_manager._cat_vectors["tissue"]._non_validated == ["pulmo"]
    assert curator._cat_manager._cat_vectors["tissue"]._synonyms == {"pulmo": "lung"}

    curator.cat.standardize("tissue")
    assert curator._cat_manager._cat_vectors["tissue"]._non_validated == []
