import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from lamindb.core.exceptions import InvalidArgument
from lamindb.curators import save_artifact


def test_nullable():
    disease = ln.Feature(name="disease", dtype=ln.ULabel, nullable=False).save()
    schema = ln.Schema(features=[disease]).save()
    dataset = {"disease": pd.Categorical([pd.NA, "asthma"])}
    df = pd.DataFrame(dataset)
    curator = ln.curators.DataFrameCurator(df, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as e:
        assert str(e).startswith("non-nullable series 'disease' contains null values")
    # make feature nullable
    # (needs to throw an error if already datasets were validated with it)
    disease.nullable = True
    disease.save()
    curator = ln.curators.DataFrameCurator(df, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as e:
        assert str(e).startswith("1 term is not validated: 'asthma'")
    schema.delete()
    disease.delete()


def test_save_artifact_invalid_data_type():
    data = np.array([1, 2, 3])

    # Check that the correct exception is raised with the expected message
    with pytest.raises(
        InvalidArgument,
        match="data must be one of pd.Dataframe, AnnData, MuData, SpatialData.",
    ):
        save_artifact(data=data, fields={"field1": "attr1"})
