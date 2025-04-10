import re

import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from lamindb.core.exceptions import InvalidArgument, ValidationError
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
    with pytest.raises(
        ValidationError,
        # match=re.escape("1 term is not validated: 'asthma'"),  # TODO: need the message
    ):
        curator.validate()
    schema.delete()
    disease.delete()


def test_save_artifact_invalid_data_type():
    data = np.array([1, 2, 3])

    with pytest.raises(
        InvalidArgument,
        match=re.escape(
            "data must be one of pd.Dataframe, AnnData, MuData, SpatialData."
        ),
    ):
        save_artifact(data=data, fields={"field1": "attr1"})


def test_pandera_dataframe_schema():
    # DataFrames
    df = pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
        }
    )
    # missing a column
    df_missing_column = pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
        }
    )
    # changed columns order
    df_changed_order = pd.DataFrame(
        {
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
            "sample_id": ["sample1", "sample2"],
        }
    )
    # additional column
    df_extra_column = pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
            "extra_column": ["Extra 1", "Extra 2"],
        }
    )

    # schemas
    schema_minimal_set = ln.Schema(
        name="my-schema minimal_set",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save(),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
    ).save()
    schema_maximal_set = ln.Schema(
        name="my-schema maximal_set",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save(),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
        minimal_set=False,
        maximal_set=True,
    ).save()
    schema_ordered_set = ln.Schema(
        name="my-schema ordered_set",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save(),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
        ordered_set=True,
    ).save()

    # minimal_set=True
    ln.curators.DataFrameCurator(df, schema=schema_minimal_set).validate()
    # can't miss a required column
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_missing_column, schema=schema_minimal_set
        ).validate()
    # doesn't care about order
    ln.curators.DataFrameCurator(df_changed_order, schema=schema_minimal_set).validate()
    # extra column is fine
    ln.curators.DataFrameCurator(df_extra_column, schema=schema_minimal_set).validate()

    # minimal_set=False
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_extra_column, schema=schema_maximal_set
        ).validate()
    ln.curators.DataFrameCurator(
        df_missing_column, schema=schema_maximal_set
    ).validate()

    # ordered_set=True
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_changed_order, schema=schema_ordered_set
        ).validate()

    # default=True (require) for a single feature when minimal_set=False
    # note this modifies the "sample_type" feature to be required
    schema_not_minimal_set_require_sample_type = ln.Schema(
        name="my-schema minimal_set require sample_type",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save(),
            ln.Feature(name="sample_type", dtype=str, required=True).save(),
        ],
        minimal_set=False,
    ).save()
    # missing "sample_type" column now raises an error
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_missing_column, schema=schema_not_minimal_set_require_sample_type
        ).validate()
    # missing a non-required column is fine
    ln.curators.DataFrameCurator(
        pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "sample_type": ["Type A", "Type B"],
            }
        ),
        schema=schema_not_minimal_set_require_sample_type,
    ).validate()

    # clean up
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
