import re

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.core.exceptions import ValidationError


@pytest.fixture
def df():
    # DataFrames
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
        }
    )


@pytest.fixture
def df_missing_sample_type_column():
    # missing a column
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
        }
    )


@pytest.fixture
def df_missing_sample_name_column():
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_type": ["Type A", "Type B"],
        }
    )


@pytest.fixture
def df_changed_order():
    # changed columns order
    return pd.DataFrame(
        {
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
            "sample_id": ["sample1", "sample2"],
        }
    )


@pytest.fixture
def df_extra_column():
    # additional column
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
            "extra_column": ["Extra 1", "Extra 2"],
        }
    )


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


def test_pandera_dataframe_schema(
    df,
    df_missing_sample_type_column,
    df_changed_order,
    df_extra_column,
    df_missing_sample_name_column,
):
    # schemas
    schema_all_required = ln.Schema(
        name="my-schema all required",
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

    # minimal_set=True, all three columns are required
    ln.curators.DataFrameCurator(df, schema=schema_all_required).validate()
    # can't miss a required column
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_missing_sample_type_column, schema=schema_all_required
        ).validate()
    # doesn't care about order
    ln.curators.DataFrameCurator(
        df_changed_order, schema=schema_all_required
    ).validate()
    # extra column is fine
    ln.curators.DataFrameCurator(df_extra_column, schema=schema_all_required).validate()

    # maximal_set=True, extra column is not allowed
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_extra_column,
            schema=schema_maximal_set,  # extra column is not allowed
        ).validate()
    # minimal_set=False, missing column is allowed
    ln.curators.DataFrameCurator(
        df_missing_sample_type_column, schema=schema_maximal_set
    ).validate()

    # ordered_set=True, order matters
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_changed_order, schema=schema_ordered_set
        ).validate()

    # a feature is optional
    schema_optional_sample_name = ln.Schema(
        name="my-schema optional sample_name",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save().with_config(optional=True),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
    ).save()
    # missing required "sample_type" column raises an error
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            df_missing_sample_type_column,
            schema=schema_optional_sample_name,
        ).validate()
    # missing optional column "sample_name" is fine
    ln.curators.DataFrameCurator(
        df_missing_sample_name_column, schema=schema_optional_sample_name
    ).validate()

    # clean up
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()


def test_schema_optionals():
    schema = ln.Schema(
        name="my-schema",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save().with_config(optional=True),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
    ).save()
    assert schema.optionals.get().list("name") == ["sample_name"]

    # set sample_type to optional
    with pytest.raises(
        TypeError,
        match=re.escape("features must be a list of Feature records!"),
    ):
        schema.optionals.set("test")
    schema.optionals.set([ln.Feature.get(name="sample_type")])
    assert schema.optionals.get().list("name") == ["sample_type"]
    # add sample_name to optionals
    with pytest.raises(
        TypeError,
        match=re.escape("features must be a list of Feature records!"),
    ):
        schema.optionals.add("test")
    schema.optionals.add(ln.Feature.get(name="sample_name"))
    assert schema.optionals.get().list("name") == ["sample_name", "sample_type"]

    # clean up
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()


def test_schema_ordered_set(df):
    # create features with a different order so that sample_id is not the first
    ln.Feature(name="sample_name", dtype=str).save()
    ln.Feature(name="sample_type", dtype=str).save()
    ln.Feature(name="sample_id", dtype=str).save()

    # create an ordered schema with sample_id as the first feature
    schema = ln.Schema(
        name="my-schema",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save(),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
        ordered_set=True,
    ).save()

    assert ln.curators.DataFrameCurator(df, schema=schema).validate() is None

    # clean up
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()


def test_schema_minimal_set_var(df):
    adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")

    var_schema = ln.Schema(
        itype=bt.Gene.ensembl_gene_id,
    ).save()
    schema = ln.Schema(otype="AnnData", slots={"var": var_schema}).save()
    curator = ln.curators.AnnDataCurator(adata, schema)
    curator.validate()


def test_schema_maximal_set_var(df):
    adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")

    var_schema = ln.Schema(itype=bt.Gene.ensembl_gene_id).save()
    schema = ln.Schema(otype="AnnData", slots={"var": var_schema}).save()
    curator = ln.curators.AnnDataCurator(adata, schema)
    curator.validate()
