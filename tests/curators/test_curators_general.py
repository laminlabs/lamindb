import re
import textwrap

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.core.exceptions import ValidationError


def _strip_ansi(text: str) -> str:
    """Remove ANSI escape sequences from a string."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


@pytest.fixture
def df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
        }
    )


@pytest.fixture
def df_missing_sample_type_column() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
        }
    )


@pytest.fixture
def df_missing_sample_name_column() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_type": ["Type A", "Type B"],
        }
    )


@pytest.fixture
def df_changed_col_order() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
            "sample_id": ["sample1", "sample2"],
        }
    )


@pytest.fixture
def df_extra_column() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
            "sample_name": ["Sample 1", "Sample 2"],
            "sample_type": ["Type A", "Type B"],
            "extra_column": ["Extra 1", "Extra 2"],
        }
    )


@pytest.fixture
def df_disease() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "disease": pd.Categorical(
                [
                    # Only after 2025 mondo
                    "HDAC4-related haploinsufficiency syndrome",
                    "SAMD9L-related spectrum and myeloid neoplasm risk",
                    # Already before 2025 mondo
                    "essential hypertension",
                    "essential hypertension",
                    "asthma",
                ]
            ),
        }
    )


@pytest.fixture
def disease_ontology_old() -> bt.Source:
    return bt.Disease.add_source(
        bt.Source.using("laminlabs/bionty-assets")
        .get(entity="bionty.Disease", version="2024-08-06", organism="all")
        .save()
    )


def test_curator__repr__(df):
    schema = ln.Schema(
        name="sample schema",
        features=[ln.Feature(name="sample_id", dtype="str").save()],
    ).save()
    curator = ln.curators.DataFrameCurator(df, schema)

    expected_repr = textwrap.dedent("""\
    DataFrameCurator(Schema: sample schema, unvalidated)
    """).strip()

    actual_repr = _strip_ansi(repr(curator))
    print(actual_repr)
    assert actual_repr.strip() == expected_repr.strip()

    schema.delete()


def test_nullable():
    disease = ln.Feature(name="disease", dtype=ln.ULabel, nullable=False).save()
    schema = ln.Schema(features=[disease]).save()
    dataset = {"disease": pd.Categorical([pd.NA, "asthma"])}
    df = pd.DataFrame(dataset)
    curator = ln.curators.DataFrameCurator(df, schema)
    with pytest.raises(ln.errors.ValidationError) as err:
        curator.validate()
    assert "non-nullable series 'disease' contains null values" in err.exconly()
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
    df_changed_col_order,
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
        df_changed_col_order, schema=schema_all_required
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
            df_changed_col_order, schema=schema_ordered_set
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


def test_schema_not_saved(df):
    """Attempting to validate an unsaved Schema must error."""
    feature = ln.Feature(name="cell_type", dtype="str").save()
    schema = ln.Schema(features=[feature])

    with pytest.raises(ValueError) as excinfo:
        ln.curators.DataFrameCurator(df, schema)
    assert excinfo.exconly() == (
        "ValueError: Schema must be saved before curation. Please save it using '.save()'."
    )


def test_schema_optionals():
    schema = ln.Schema(
        name="my-schema",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save().with_config(optional=True),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
    ).save()
    assert schema.optionals.get().list("name") == [
        "sample_name",
    ]

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


@pytest.mark.parametrize("minimal_set", [True, False])
def test_schema_minimal_set_var_allowed(minimal_set):
    """Independent of the value of minimal_set, invalid ensembl gene IDs are allowed."""
    adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")
    adata.var_names = [adata.var_names[0], adata.var_names[1], "NOT_VALID_ENSEMBL"]

    var_schema = ln.Schema(
        itype=bt.Gene.ensembl_gene_id,
        minimal_set=minimal_set,
    ).save()
    schema = ln.Schema(otype="AnnData", slots={"var.T": var_schema}).save()
    curator = ln.curators.AnnDataCurator(adata, schema)
    curator.validate()

    # clean up
    schema.delete()


def test_schema_maximal_set_var():
    """If maximal_set is True, invalid ensembl gene IDs are not allowed."""
    adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")
    adata.var_names = [adata.var_names[0], adata.var_names[1], "NOT_VALID_ENSEMBL"]

    var_schema = ln.Schema(itype=bt.Gene.ensembl_gene_id, maximal_set=True).save()
    schema = ln.Schema(otype="AnnData", slots={"var.T": var_schema}).save()

    curator = ln.curators.AnnDataCurator(adata, schema)
    with pytest.raises(ValidationError) as error:
        curator.validate()
    assert error.exconly() == (
        "lamindb.errors.ValidationError: 1 term not validated in feature 'columns' in slot 'var.T': 'NOT_VALID_ENSEMBL'\n"
        "    → fix typos, remove non-existent values, or save terms via: curator.slots['var.T'].cat.add_new_from('columns')"
    )

    # clean up
    schema.delete()


def test_feature_dtype_path():
    df = pd.DataFrame(
        {
            "sample": ["Sample_X", "Sample_Y", "Sample_Y"],
            "fastq_1": [
                "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R1_001.fastq.gz",
                "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L001_R1_001.fastq.gz",
                "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L002_R1_001.fastq.gz",
            ],
            "fastq_2": [
                "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R2_001.fastq.gz",
                "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L001_R2_001.fastq.gz",
                "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L002_R2_001.fastq.gz",
            ],
            "expected_cells": [5000, 5000, 5000],
        }
    )

    nextflow_schema = ln.Schema(
        name="nf-core/scrnaseq pipeline - params.input schema",
        description="https://github.com/nf-core/scrnaseq/blob/4.0.0/assets/schema_input.json",
        features=[
            ln.Feature(
                name="sample",
                dtype="str",
                nullable=False,
                description="Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (_).",
            ).save(),
            ln.Feature(
                name="fastq_1",
                dtype="path",
                nullable=False,
                description="Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension “.fastq.gz” or “.fq.gz”.",
            ).save(),
            ln.Feature(
                name="fastq_2",
                dtype="path",
                description="Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension “.fastq.gz” or “.fq.gz”.",
            ).save(),
            ln.Feature(
                name="expected_cells",
                dtype=int,
                description="Number of cells expected for a sample. Must be an integer. If multiple rows are provided for the same sample, this must be the same number for all rows, i.e. the total number of expected cells for the sample.",
            ).save(),
            ln.Feature(
                name="seq_center",
                dtype=str,
                description="Sequencing center for the sample. If multiple rows are provided for the same sample, this must be the same string for all rows. Samples sequenced at different centers are considered different samples and must have different identifiers.",
            ).save(),
            ln.Feature(
                name="sample_type",
                dtype=str,
                description='"atac", "gex"',
            ).save(),
            ln.Feature(
                name="feature_type",
                dtype=str,
                description='"gex", "vdj", "ab", "crispr", "cmo"',
            ).save(),
        ],
    ).save()

    nextflow_schema.optionals.set(
        [
            ln.Feature.get(name="expected_cells"),
            ln.Feature.get(name="seq_center"),
            ln.Feature.get(name="sample_type"),
            ln.Feature.get(name="feature_type"),
        ]
    )

    curator = ln.curators.DataFrameCurator(df, schema=nextflow_schema)
    assert curator.validate() is None

    # clean up
    nextflow_schema.delete()
    ln.Feature.filter().delete()


def test_cat_filters_specific_source_uid(df_disease, disease_ontology_old):
    """Specific source_uid passed to the `cat_filters`"""
    schema = ln.Schema(
        features=[
            ln.Feature(
                name="disease",
                dtype=bt.Disease,
                cat_filters={"source__uid": disease_ontology_old.uid},
            ).save(),
        ],
    ).save()

    curator = ln.curators.DataFrameCurator(df_disease, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert (
            "2 terms not validated in feature 'disease': 'HDAC4-related haploinsufficiency syndrome', 'SAMD9L-related spectrum and myeloid neoplasm risk'"
            in str(error)
        )

    schema.delete()


def test_cat_filters_specific_source(df_disease, disease_ontology_old):
    """Specific Source record passed to the `cat_filters`"""
    schema = ln.Schema(
        features=[
            ln.Feature(
                name="disease",
                dtype=bt.Disease,
                cat_filters={"source": disease_ontology_old},
            ).save(),
        ],
    ).save()

    curator = ln.curators.DataFrameCurator(df_disease, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert (
            "2 terms not validated in feature 'disease': 'HDAC4-related haploinsufficiency syndrome', 'SAMD9L-related spectrum and myeloid neoplasm risk'"
            in str(error)
        )

    schema.delete()


def test_cat_filters_multiple_relation_filters(df_disease, disease_ontology_old):
    """Multiple relation filters in cat_filters"""
    schema = ln.Schema(
        features=[
            ln.Feature(
                name="disease",
                dtype=bt.Disease,
                cat_filters={
                    "source__uid": disease_ontology_old.uid,
                    "organism__name": "all",
                },
            ).save(),
        ],
    ).save()
    curator = ln.curators.DataFrameCurator(df_disease, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert (
            "2 terms not validated in feature 'disease': 'HDAC4-related haploinsufficiency syndrome', 'SAMD9L-related spectrum and myeloid neoplasm risk'"
            in str(error)
        )
    schema.delete()


def test_curate_columns(df):
    """Test that columns can be curated."""
    schema = ln.Schema(
        name="sample schema",
        features=[
            ln.Feature(name="sample_id", dtype="str").save(),
            ln.Feature(name="sample_name", dtype="str").save(),
            ln.Feature(name="sample_type", dtype="str").save(),
        ],
    ).save()

    # make one column name invalid
    df.rename(columns={"sample_name": "sample_name_name"}, inplace=True)

    curator = ln.curators.DataFrameCurator(df, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert "column 'sample_name' not in dataframe" in str(error)

    # now fix the column
    df.rename(columns={"sample_name_name": "sample_name"}, inplace=True)
    curator.validate()

    schema.delete()
    ln.Feature.filter().delete()


def test_wrong_datatype(df):
    feature = ln.Feature(name="sample_id", dtype=ln.ULabel).save()
    schema = ln.Schema(features=[feature]).save()

    curator = ln.curators.DataFrameCurator(df, schema)
    with pytest.raises(ln.errors.ValidationError) as excinfo:
        curator.validate()

    assert "expected series 'sample_id' to have type category, got object" in str(
        excinfo.value
    )
    assert (
        "Hint: Consider setting 'coerce_datatype=True' to attempt coercing/converting values during validation to the pre-defined dtype."
        in str(excinfo.value)
    )

    schema.delete()
    feature.delete()


def test_hash_index_feature(df):
    df_index = df.set_index("sample_id")
    sample_name = ln.Feature(name="sample_name", dtype="str").save()
    sample_name.uid = "OpQAD5Ifu89t"
    sample_name.save()
    sample_type = ln.Feature(name="sample_type", dtype="str").save()
    sample_type.uid = "7I4u69RiCAVy"
    sample_type.save()
    sample_id = ln.Feature(name="sample_id", dtype="str").save()
    sample_id.uid = "uValv1YfEQib"
    sample_id.save()
    schema_index = ln.Schema(
        name="sample schema with index",
        features=[
            sample_name,
            sample_type,
        ],
        index=sample_id,
    ).save()
    assert schema_index.hash == "drtQMP4N4xEebS49DO-9Jw"

    schema = ln.Schema(
        name="sample schema",
        features=[
            sample_id,
            sample_name,
            sample_type,
        ],
    ).save()
    assert schema.hash == "Z_dmk1WendD15s2FyBW1HA"

    artifact = ln.Artifact.from_df(
        df_index, key="curated_df.parquet", schema=schema_index
    ).save()
    assert artifact.feature_sets.all().one() == schema_index

    # clean up
    artifact.delete(permanent=True)
    schema_index.delete()
    schema.delete()
    ln.Feature.filter().delete()


def test_add_new_from_subtype(df):
    """Test that add_new_from works with subtypes."""
    sample_type = ln.Record(name="SampleType", is_type=True).save()
    ln.Record(name="Type A", type=sample_type).save()
    schema = ln.Schema(
        name="sample schema",
        features=[
            ln.Feature(name="sample_id", dtype="str").save(),
            ln.Feature(name="sample_name", dtype="str").save(),
            ln.Feature(name="sample_type", dtype="cat[Record[SampleType]]").save(),
        ],
        coerce_dtype=True,
    ).save()

    curator = ln.curators.DataFrameCurator(df, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert "1 term not validated in feature 'sample_type': 'Type B'" in str(error)

    # add new from subtype
    curator.cat.add_new_from("sample_type")
    curator.validate()
    assert sample_type.records.list("name") == ["Type A", "Type B"]

    # clean up
    schema.delete()
    ln.Feature.filter().delete()
    ln.Record.filter().update(type=None)
    ln.Record.filter().delete()
