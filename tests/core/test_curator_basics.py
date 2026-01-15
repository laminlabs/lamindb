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
        bt.Source.connect("laminlabs/bionty-assets")
        .get(entity="bionty.Disease", version="2024-08-06", organism="all")
        .save()
    )


@pytest.fixture(scope="module")
def lists_df():
    return pd.DataFrame(
        {
            "sample_id": [["sample1", "sample2"], ["sample2"], ["sample3"]],
            "dose": [[1.2, 2.3], [1.2], [2.3]],
            "cell_type": [["B cell", "T cell"], ["B cell"], ["T cell"]],
            "tissue": [["blood", "pulmo"], ["blood"], ["lung"]],
        }
    )


@pytest.fixture(scope="module")
def cat_df():
    return pd.DataFrame(
        {
            "sample_id": [["sample1", "sample2"], ["sample2"], ["sample3"]],
            "dose": [[1.2, 2.3], [1.2], [2.3]],
            "cell_type": [["B cell", "T cell"], ["B cell"], ["T cell"]],
            "tissue": ["blood", "blood", "lung"],
        }
    )


def test_curator_df_multivalue(lists_df, cat_df):
    feature1 = ln.Feature(name="sample_id", dtype=list[str]).save()
    feature2 = ln.Feature(name="dose", dtype=list[float]).save()
    feature3 = ln.Feature(name="cell_type", dtype=list[str]).save()
    feature4 = ln.Feature(name="tissue", dtype=list[bt.Tissue]).save()
    schema = ln.Schema(
        name="lists schema cat",
        features=[
            feature1,
            feature2,
            feature3,
            feature4,
        ],
    ).save()

    curator = ln.curators.DataFrameCurator(lists_df, schema)
    with pytest.raises(ValidationError):
        curator.validate()
    assert curator.cat._cat_vectors.keys() == {"columns", "tissue"}
    assert curator.cat._cat_vectors["tissue"]._validated == ["blood", "lung"]
    assert curator.cat._cat_vectors["tissue"]._non_validated == ["pulmo"]
    assert curator.cat._cat_vectors["tissue"]._synonyms == {"pulmo": "lung"}

    curator.cat.standardize("tissue")
    assert curator.cat._cat_vectors["tissue"]._non_validated == []
    assert lists_df["tissue"].tolist() == [["blood", "lung"], ["blood"], ["lung"]]

    assert curator.validate() is None

    # test with cat_df which has a non-list tissue
    curator = ln.curators.DataFrameCurator(cat_df, schema)
    with pytest.raises(ValidationError):
        curator.validate()

    schema.delete(permanent=True)
    feature1.delete(permanent=True)
    feature2.delete(permanent=True)
    feature3.delete(permanent=True)
    feature4.delete(permanent=True)


def test_curators_list_feature_nullable_empty_list():
    """Test that a list feature that is nullable can accept empty lists."""
    feature_list = ln.Feature(
        name="list_tissue", dtype=list[bt.Tissue.ontology_id], nullable=True
    ).save()
    feature_int = ln.Feature(name="feature int", dtype=int, nullable=True).save()
    schema = ln.Schema(
        name="test_list_feature_schema",
        features=[feature_list, feature_int],
        coerce=True,
    ).save()

    df = pd.DataFrame({"list_tissue": [], "feature int": []})
    ln.curators.DataFrameCurator(df, schema).validate()

    # clean up
    schema.delete(permanent=True)
    feature_list.delete(permanent=True)
    feature_int.delete(permanent=True)


def test_curator__repr__(df):
    feature = ln.Feature(name="sample_id", dtype="str").save()
    schema = ln.Schema(
        name="sample schema",
        features=[feature],
    ).save()
    curator = ln.curators.DataFrameCurator(df, schema)

    expected_repr = textwrap.dedent("""\
    DataFrameCurator(Schema: sample schema, unvalidated)
    """).strip()

    actual_repr = _strip_ansi(repr(curator))
    print(actual_repr)
    assert actual_repr.strip() == expected_repr.strip()

    schema.delete(permanent=True)
    feature.delete(permanent=True)


@pytest.mark.parametrize(
    "model_class",
    [ln.ULabel, ln.Record],
)
def test_df_curator_typed_categorical(model_class):
    # root level
    sample_root_type = model_class(name="Sample", is_type=True).save()
    for name in ["s1", "s2"]:
        model_class(name=name, type=sample_root_type).save()

    # lab A level
    lab_a_type = model_class(name="LabA", is_type=True).save()
    sample_a_type = model_class(name="Sample", is_type=True, type=lab_a_type).save()
    for name in ["s3", "s4"]:
        model_class(name=name, type=sample_a_type).save()

    # lab B level
    lab_b_type = model_class(name="LabB", is_type=True).save()
    sample_b_type = model_class(name="Sample", is_type=True, type=lab_b_type).save()
    for name in ["s5", "s6"]:
        model_class(name=name, type=sample_b_type).save()

    df = pd.DataFrame(
        {
            "biosample_name": pd.Categorical(["s1", "s2", "s3", "s4", "s5", "s6"]),
        }
    )

    feature = ln.Feature(name="biosample_name", dtype=sample_a_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert "4 terms not validated in feature 'biosample_name':" in error.exconly()
    assert set(curator.cat._cat_vectors["biosample_name"]._validated) == {
        "s3",
        "s4",
    }
    assert set(curator.cat._cat_vectors["biosample_name"]._non_validated) == {
        "s1",
        "s2",
        "s5",
        "s6",
    }

    # Move LabB under LabA
    lab_b_type.type = lab_a_type
    lab_b_type.save()
    feature.delete(permanent=True)  # re-create the feature with the new dtype
    feature = ln.Feature(name="biosample_name", dtype=lab_a_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert set(curator.cat._cat_vectors["biosample_name"]._validated) == {
        "s3",
        "s4",
        "s5",
        "s6",
    }
    assert set(curator.cat._cat_vectors["biosample_name"]._non_validated) == {
        "s1",
        "s2",
    }

    # Lab at the root
    feature.delete(permanent=True)  # re-create the feature with the new dtype
    feature = ln.Feature(name="biosample_name", dtype=sample_root_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert set(curator.cat._cat_vectors["biosample_name"]._validated) == {
        "s1",
        "s2",
    }
    assert set(curator.cat._cat_vectors["biosample_name"]._non_validated) == {
        "s3",
        "s4",
        "s5",
        "s6",
    }

    attribute = model_class.__name__.lower() + "s"
    getattr(sample_a_type, attribute).all().delete(permanent=True)
    getattr(sample_b_type, attribute).all().delete(permanent=True)
    getattr(lab_b_type, attribute).all().delete(permanent=True)
    getattr(lab_a_type, attribute).all().delete(permanent=True)
    lab_a_type.delete(permanent=True)
    lab_b_type.delete(permanent=True)
    getattr(sample_root_type, attribute).all().delete(permanent=True)
    sample_root_type.delete(permanent=True)
    feature.delete(permanent=True)


def test_df_curator_same_name_at_different_levels_involving_root():
    s1_root = ln.Record(name="s1").save()
    lab_a_type = ln.Record(name="LabA", is_type=True).save()
    s1_lab_a = ln.Record(name="s1", type=lab_a_type).save()
    df = pd.DataFrame({"biosample_name": pd.Categorical(["s1"])})

    # feature constraining to lab_a_type
    feature = ln.Feature(name="biosample_name", dtype=lab_a_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    curator.validate()
    cat_vector = curator._atomic_curator.cat._cat_vectors["biosample_name"]
    assert cat_vector._validated == ["s1"]
    assert len(cat_vector.records) == 1
    assert cat_vector.records[0] == s1_lab_a

    # feature constraining to root
    feature.delete(permanent=True)
    feature = ln.Feature(name="biosample_name", dtype=ln.Record).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    curator.validate()
    cat_vector = curator._atomic_curator.cat._cat_vectors["biosample_name"]
    assert cat_vector._validated == ["s1"]
    assert len(cat_vector.records) == 1
    assert cat_vector.records[0] == s1_root

    feature.delete(permanent=True)
    s1_root.delete(permanent=True)
    s1_lab_a.delete(permanent=True)
    lab_a_type.delete(permanent=True)


def test_df_curator_same_name_at_different_levels_below_root():
    department_a_type = ln.Record(name="DepartmentA", is_type=True).save()
    s1_department_a = ln.Record(name="s1", type=department_a_type).save()
    lab_a_type = ln.Record(name="LabA", is_type=True, type=department_a_type).save()
    s1_lab_a = ln.Record(name="s1", type=lab_a_type).save()
    df = pd.DataFrame({"biosample_name": pd.Categorical(["s1"])})

    # feature constraining to lab_a_type
    feature = ln.Feature(name="biosample_name", dtype=lab_a_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    curator.validate()
    cat_vector = curator._atomic_curator.cat._cat_vectors["biosample_name"]
    assert cat_vector._validated == ["s1"]
    assert len(cat_vector.records) == 1
    assert cat_vector.records[0] == s1_lab_a

    # feature constraining to department_a_type
    feature.delete(permanent=True)
    feature = ln.Feature(name="biosample_name", dtype=department_a_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    curator.validate()
    cat_vector = curator._atomic_curator.cat._cat_vectors["biosample_name"]
    assert cat_vector._validated == ["s1"]
    assert len(cat_vector.records) == 1
    assert cat_vector.records[0] == s1_department_a

    feature.delete(permanent=True)
    s1_department_a.delete(permanent=True)
    s1_lab_a.delete(permanent=True)
    lab_a_type.delete(permanent=True)
    department_a_type.delete(permanent=True)


def test_df_curator_same_name_at_same_level():
    # below root level
    lab_a_type = ln.Record(name="LabA", is_type=True).save()
    record_1 = ln.Record(name="s1", type=lab_a_type).save()
    lab_b_type = ln.Record(name="LabB", is_type=True).save()
    record_2 = ln.Record(name="s1", type=lab_b_type).save()
    df = pd.DataFrame({"biosample_name": pd.Categorical(["s1"])})
    feature = ln.Feature(name="biosample_name", dtype=ln.Record).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        "Ambiguous match for Record 's1': found 2 records at depth 1 (under types: ['LabA', 'LabB'])"
        in error.exconly()
    )

    # at root level
    record_1.type = None
    record_1.save()
    record_2.type = None
    record_2.save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        "Ambiguous match for Record 's1': found 2 root-level records" in error.exconly()
    )

    feature.delete(permanent=True)
    record_1.delete(permanent=True)
    lab_a_type.delete(permanent=True)
    record_2.delete(permanent=True)
    lab_b_type.delete(permanent=True)


# also see test_features_name_duplicates_across_equal_levels
def test_curator_schema_feature_mapping():
    lab_a_type = ln.Feature(name="LabA", is_type=True).save()
    feature1 = ln.Feature(name="sample_name", dtype="str", type=lab_a_type).save()
    lab_b_type = ln.Feature(name="LabB", is_type=True).save()
    feature2 = ln.Feature(name="sample_name", dtype="str", type=lab_b_type).save()
    schema = ln.Schema([feature1], name="Lab A schema").save()
    df = pd.DataFrame({"sample_name": ["s1", "s2"]})
    curator = ln.curators.DataFrameCurator(df, schema)
    curator.validate()
    cat_vector = curator._atomic_curator.cat._cat_vectors["columns"]
    assert len(cat_vector.records) == 1
    assert len(cat_vector._validated) == 1
    schema.delete(permanent=True)
    feature1.delete(permanent=True)
    feature2.delete(permanent=True)
    lab_a_type.delete(permanent=True)
    lab_b_type.delete(permanent=True)


def test_dtypes_at_different_levels(ccaplog):
    sample_type_root = ln.Record(name="Sample", is_type=True).save()
    lab_a_type = ln.Record(name="LabA", is_type=True).save()
    sample_type_a = ln.Record(name="Sample", is_type=True, type=lab_a_type).save()
    s1_lab_a = ln.Record(name="s1", type=sample_type_a).save()
    df = pd.DataFrame({"biosample_name": pd.Categorical(["s1"])})
    feature = ln.Feature(name="biosample_name", dtype=sample_type_root).save()
    schema = ln.Schema(features=[feature]).save()
    sample_type_root.delete()
    df = pd.DataFrame({"biosample_name": pd.Categorical(["s1"])})
    # UID-based lookup can find records in trash, so curator creation should succeed
    # but a warning should be printed
    curator = ln.curators.DataFrameCurator(df, schema)
    assert "from trash" in ccaplog.text
    schema.delete(permanent=True)
    sample_type_root.restore()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert "1 term not validated in feature 'biosample_name': 's1'" in error.exconly()
    s1_root = ln.Record(name="s1", type=sample_type_root).save()
    curator.validate()
    cat_vector = curator._atomic_curator.cat._cat_vectors["biosample_name"]
    assert cat_vector._validated == ["s1"]
    assert len(cat_vector.records) == 1
    assert cat_vector.records[0] == s1_root
    # update feature dtype
    feature.delete(permanent=True)
    feature = ln.Feature(name="biosample_name", dtype=sample_type_a).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    curator.validate()
    cat_vector = curator._atomic_curator.cat._cat_vectors["biosample_name"]
    assert cat_vector._validated == ["s1"]
    assert len(cat_vector.records) == 1
    assert cat_vector.records[0] == s1_lab_a
    feature.delete(permanent=True)
    s1_lab_a.delete(permanent=True)
    sample_type_a.delete(permanent=True)
    lab_a_type.delete(permanent=True)
    s1_root.delete(permanent=True)
    sample_type_root.delete(permanent=True)


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

    schema.delete(permanent=True)
    disease.delete(permanent=True)


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

    # maximal_set=True, extra column is *not* allowed
    # check that __lamindb values are OK
    df["__lamindb_record_uid__"] = "some_value"
    ln.curators.DataFrameCurator(df, schema=schema_maximal_set).validate()
    del df["__lamindb_record_uid__"]
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
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


def test_schema_not_saved(df):
    """Attempting to validate an unsaved Schema must error."""
    feature = ln.Feature(name="cell_type", dtype=str).save()
    schema = ln.Schema(features=[feature])

    with pytest.raises(ValueError) as excinfo:
        ln.curators.DataFrameCurator(df, schema)
    assert excinfo.exconly() == (
        "ValueError: Schema must be saved before curation. Please save it using '.save()'."
    )


def test_schema_artifact_annotated(df):
    """A passed Artifact should be annotated with a Schema if successfully curated."""
    af = ln.Artifact.from_dataframe(df, key="test.parquet").save()
    schema = ln.Schema(
        name="sample schema",
        features=[ln.Feature(name="sample_id", dtype="str").save()],
    ).save()
    curator = ln.curators.DataFrameCurator(af, schema)
    curator.validate()
    curator.save_artifact()
    af_queried = ln.Artifact.filter(key="test.parquet").one()
    assert af_queried.schema is not None

    # clean up
    af.delete(permanent=True)
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


def test_schema_optionals():
    schema = ln.Schema(
        name="my-schema",
        features=[
            ln.Feature(name="sample_id", dtype=str).save(),
            ln.Feature(name="sample_name", dtype=str).save().with_config(optional=True),
            ln.Feature(name="sample_type", dtype=str).save(),
        ],
    ).save()
    assert schema.optionals.get().to_list("name") == [
        "sample_name",
    ]

    # set sample_type to optional
    with pytest.raises(
        TypeError,
        match=re.escape("features must be a list of Feature records!"),
    ):
        schema.optionals.set("test")
    schema.optionals.set([ln.Feature.get(name="sample_type")])
    assert schema.optionals.get().to_list("name") == ["sample_type"]
    # add sample_name to optionals
    with pytest.raises(
        TypeError,
        match=re.escape("features must be a list of Feature records!"),
    ):
        schema.optionals.add("test")
    schema.optionals.add(ln.Feature.get(name="sample_name"))
    assert schema.optionals.get().to_list("name") == ["sample_name", "sample_type"]

    # clean up
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


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
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


@pytest.mark.parametrize("minimal_set", [True, False])
def test_schema_minimal_set_var_allowed(minimal_set):
    """Independent of the value of minimal_set, invalid ensembl gene IDs are allowed."""
    adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
    adata.var_names = [adata.var_names[0], adata.var_names[1], "NOT_VALID_ENSEMBL"]

    var_schema = ln.Schema(
        itype=bt.Gene.ensembl_gene_id,
        minimal_set=minimal_set,
    ).save()
    schema = ln.Schema(otype="AnnData", slots={"var.T": var_schema}).save()
    curator = ln.curators.AnnDataCurator(adata, schema)
    curator.validate()

    # clean up
    schema.delete(permanent=True)


def test_schema_maximal_set_var():
    """If maximal_set is True, invalid ensembl gene IDs are not allowed."""
    adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
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
    schema.delete(permanent=True)


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
    nextflow_schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


def test_cat_filters_specific_source_uid(df_disease, disease_ontology_old):
    """Specific source_uid passed to the `cat_filters`"""
    feature = ln.Feature(
        name="disease",
        dtype=bt.Disease,
        cat_filters={"source__uid": disease_ontology_old.uid},
    ).save()
    schema = ln.Schema([feature], name="test schema").save()
    curator = ln.curators.DataFrameCurator(df_disease, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert (
            "2 terms not validated in feature 'disease': 'HDAC4-related haploinsufficiency syndrome', 'SAMD9L-related spectrum and myeloid neoplasm risk'"
            in str(error)
        )
    schema.delete(permanent=True)
    feature.delete(permanent=True)


def test_cat_filters_specific_source(df_disease, disease_ontology_old):
    """Specific Source record passed to the `cat_filters`"""
    feature = ln.Feature(
        name="disease",
        dtype=bt.Disease,
        cat_filters={"source": disease_ontology_old},
    ).save()
    schema = ln.Schema([feature], name="test schema").save()
    curator = ln.curators.DataFrameCurator(df_disease, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert (
            "2 terms not validated in feature 'disease': 'HDAC4-related haploinsufficiency syndrome', 'SAMD9L-related spectrum and myeloid neoplasm risk'"
            in str(error)
        )

    schema.delete(permanent=True)
    feature.delete(permanent=True)


def test_cat_filters_multiple_relation_filters(df_disease, disease_ontology_old):
    """Multiple relation filters in cat_filters"""
    # TODO: needs to also work if both filters are from the same related model!!!
    feature = ln.Feature(
        name="disease",
        dtype=bt.Disease,
        cat_filters={
            "source__uid": disease_ontology_old.uid,
            "created_by__handle": ln.setup.settings.user.handle,
        },
    ).save()
    schema = ln.Schema([feature], name="test schema").save()
    curator = ln.curators.DataFrameCurator(df_disease, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert (
            "2 terms not validated in feature 'disease': 'HDAC4-related haploinsufficiency syndrome', 'SAMD9L-related spectrum and myeloid neoplasm risk'"
            in str(error)
        )
    schema.delete(permanent=True)
    feature.delete(permanent=True)


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

    schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


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
        "Hint: Consider setting `feature.coerce = True` to attempt coercing values during validation to the required dtype."
        in str(excinfo.value)
    )

    schema.delete(permanent=True)
    feature.delete(permanent=True)


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

    artifact = ln.Artifact.from_dataframe(
        df_index, key="curated_df.parquet", schema=schema_index
    ).save()
    assert artifact.schemas.all().one() == schema_index

    # clean up
    artifact.delete(permanent=True)
    schema_index.delete(permanent=True)
    schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


def test_add_new_from_subtype(df):
    """Test that add_new_from works with subtypes."""
    sample_type = ln.Record(name="SampleType", is_type=True).save()
    ln.Record(name="Type A", type=sample_type).save()
    schema = ln.Schema(
        name="sample schema",
        features=[
            ln.Feature(name="sample_id", dtype="str").save(),
            ln.Feature(name="sample_name", dtype="str").save(),
            ln.Feature(name="sample_type", dtype=sample_type).save(),
        ],
        coerce=True,
    ).save()

    curator = ln.curators.DataFrameCurator(df, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert "1 term not validated in feature 'sample_type': 'Type B'" in str(error)

    # add new from subtype
    curator.cat.non_validated["sample_type"]
    curator.cat.add_new_from("sample_type")
    curator.validate()
    assert sample_type.records.to_list("name") == ["Type A", "Type B"]

    # clean up
    schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    ln.Record.filter().update(type=None)
    ln.Record.filter().delete(permanent=True)


def test_index_feature_exclusion_from_categoricals(df):
    df_indexed = df.set_index("sample_id")

    sample_type_feature = ln.Feature(name="sample_type", dtype="cat[ULabel]").save()
    sample_id_feature = ln.Feature(name="sample_id", dtype="cat[ULabel]").save()

    # schema with sample_id as index (not in features)
    schema = ln.Schema(features=[sample_type_feature], index=sample_id_feature).save()

    curator = ln.curators.DataFrameCurator(df_indexed, schema)

    # Verify that only sample_type is in categoricals, not sample_id (index)
    categoricals_names = [
        f.name for f in curator._atomic_curator._cat_manager._categoricals
    ]
    assert "sample_type" in categoricals_names
    assert "sample_id" not in categoricals_names

    # Verify the cat_vectors do not include the index feature
    cat_vector_keys = list(curator.cat._cat_vectors.keys())
    assert "sample_type" in cat_vector_keys
    assert "sample_id" not in cat_vector_keys
    assert "columns" in cat_vector_keys

    # clean up
    schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
