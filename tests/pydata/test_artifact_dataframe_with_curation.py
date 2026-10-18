# Note: Almost all logic for schema-based validation is handled in the curators test suite
# This here only covers external feature annotation and validation

from pathlib import Path

import lamindb as ln
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def two_internal_features():
    feat1 = ln.Feature(name="feat1", dtype=int).save()
    feat2 = ln.Feature(name="feat2", dtype=int).save()
    yield feat1, feat2
    feat1.delete(permanent=True)
    feat2.delete(permanent=True)


@pytest.fixture(scope="module")
def two_external_features():
    feature_a = ln.Feature(name="feature_a", dtype=str).save()
    feature_b = ln.Feature(name="feature_b", dtype=str).save()
    yield feature_a, feature_b
    feature_a.delete(permanent=True)
    feature_b.delete(permanent=True)


@pytest.mark.parametrize("use_schema", [True, False])
def test_create_artifact_with_external_feature_annotations(
    use_schema: bool,
    two_external_features: tuple[ln.Feature, ln.Feature],
):
    feat1, feat2 = two_external_features
    if use_schema:
        schema = ln.Schema(features=[feat1, feat2]).save()
    else:
        schema = None
    artifact = ln.Artifact(
        ".gitignore",
        key="test_file",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    assert artifact.schema == schema
    # repeat to check idempotency (requires set_values() instead of add_values())
    artifact = ln.Artifact(
        ".gitignore",
        key="test_file",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    assert artifact.schema == schema
    if use_schema:
        with pytest.raises(ValueError) as error:
            artifact.features.remove_values("feature_a", value="x")
        assert (
            "Cannot remove values if artifact has external schema." in error.exconly()
        )
    else:
        artifact.features.remove_values("feature_a", value="x")
        assert artifact.features.get_values() == {"feature_b": "y"}
    artifact.delete(permanent=True)
    if use_schema:
        schema.delete(permanent=True)


def test_artifact_from_dataframe_with_schema(example_dataframe: pd.DataFrame):
    df = example_dataframe
    feat1 = ln.Feature(name="feat1", dtype=int).save()
    artifact = ln.Artifact.from_dataframe(
        df, key="test_df.parquet", schema="valid_features"
    ).save()
    # repeat to check idempotency
    artifact = ln.Artifact.from_dataframe(
        df, key="test_df.parquet", schema="valid_features"
    ).save()
    assert artifact.schema == ln.examples.schemas.valid_features()
    assert artifact.features.get_values() == {}
    assert (
        artifact.features.describe(return_str=True)
        == """\
Artifact: test_df.parquet (0000)
└── Dataset features
    └── columns (1)
        feat1               int"""
    )
    inferred_schema_link = artifact.schemas.through.get(artifact_id=artifact.id)
    assert inferred_schema_link.slot == "columns"
    assert inferred_schema_link.schema.members.count() == 1
    assert inferred_schema_link.schema.members.first() == feat1
    inferred_schema = inferred_schema_link.schema
    inferred_schema_link.delete()
    inferred_schema.delete(permanent=True)
    feat1.delete(permanent=True)
    artifact.delete(permanent=True)


def test_artifact_dataframe_with_features(example_dataframe: pd.DataFrame):
    """Test column names encoding when features with the same names are present."""
    artifact = ln.Artifact.from_dataframe(example_dataframe, key="df.parquet").save()
    id_feature = ln.Feature(name="id", dtype=int).save()
    uid_feature = ln.Feature(name="uid", dtype=str).save()
    artifact.features.add_values({"id": 1, "uid": "test-uid"})
    df = ln.Artifact.filter(key="df.parquet").to_dataframe(
        include=["description"], features=True
    )
    assert df.index.name == "__lamindb_artifact_id__"
    assert df.columns.tolist() == [
        "__lamindb_artifact_uid__",
        "key",
        "id",
        "uid",
        "description",
    ]
    assert df.iloc[0]["id"] == 1
    assert df.iloc[0]["uid"] == "test-uid"

    artifact.delete(permanent=True)
    id_feature.delete(permanent=True)
    uid_feature.delete(permanent=True)


def test_from_dataframe_with_external_schema(
    example_dataframe: pd.DataFrame,
    two_external_features: tuple[ln.Feature, ln.Feature],
    two_internal_features: tuple[ln.Feature, ln.Feature],
):
    df = example_dataframe
    feat1, feat2 = two_internal_features
    featA, featB = two_external_features
    schema_external = ln.Schema(features=[featA, featB]).save()

    # Case 1: wrong internal features for this dataframe
    schema_with_mistake = ln.Schema(
        features=[featA, featB],
        slots={"__external__": schema_external},
        otype="DataFrame",
    ).save()
    with pytest.raises(ln.errors.ValidationError) as error:
        artifact = ln.Artifact.from_dataframe(
            df,
            key="test_df_with_external_features.parquet",
            features={"feature_a": "x", "feature_b": "y"},
            schema=schema_with_mistake,
        ).save()
    assert "COLUMN_NOT_IN_DATAFRAME" in error.exconly()

    # alternative via DataFrameCurator directly
    with pytest.raises(ln.errors.ValidationError) as error:
        ln.curators.DataFrameCurator(
            df,
            schema=schema_with_mistake,
        ).validate()
    assert "COLUMN_NOT_IN_DATAFRAME" in error.exconly()

    # Case 2: no schema for external features provided
    schema_no_external = ln.Schema(features=[feat1, feat2]).save()
    artifact = ln.Artifact.from_dataframe(
        df,
        key="test_df_with_external_features.parquet",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema_no_external,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    artifact.delete(permanent=True)

    # alternative via DataFrameCurator directly
    curator = ln.curators.DataFrameCurator(
        df,
        schema=schema_no_external,
        features={"feature_a": "x", "feature_b": "y"},
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    artifact.delete(permanent=True)

    # Case 3: correct external schema
    schema_correct_external = ln.Schema(
        features=[feat1, feat2],
        slots={"__external__": schema_external},
        otype="DataFrame",
    ).save()

    # Case 3a: user passes no external features
    with pytest.raises(ln.errors.ValidationError) as error:
        artifact = ln.Artifact.from_dataframe(
            df,
            key="test_df_with_external_features.parquet",
            schema=schema_correct_external,
        ).save()
    assert (
        "External features slot is defined in schema but no external features were provided."
        in error.exconly()
    )

    # alternative via DataFrameCurator directly
    with pytest.raises(ln.errors.ValidationError) as error:
        curator = ln.curators.DataFrameCurator(
            df,
            schema=schema_correct_external,
        )
        artifact = curator.save_artifact(
            key="test_df_with_external_features.parquet",
        ).save()
    assert (
        "External features slot is defined in schema but no external features were provided."
        in error.exconly()
    )

    # Case 3b: user provides external features
    artifact = ln.Artifact.from_dataframe(
        df,
        key="test_df_with_external_features.parquet",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema_correct_external,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    assert (
        artifact.features.describe(return_str=True)
        == """\
Artifact: test_df_with_external_features.parquet (0000)
├── Dataset features
│   └── columns (2)
│       feat1               int
│       feat2               int
└── External features
    └── feature_a           str                      x
        feature_b           str                      y"""
    )
    with pytest.raises(ValueError) as error:
        artifact.features.remove_values("feature_a", value="x")
    assert "Cannot remove values if artifact has external schema." in error.exconly()
    artifact.delete(permanent=True)

    # alternative via DataFrameCurator directly
    curator = ln.curators.DataFrameCurator(
        df,
        schema=schema_correct_external,
        features={"feature_a": "x", "feature_b": "y"},
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}

    # call this again to check calling with an existing artifact
    curator = ln.curators.DataFrameCurator(
        artifact,
        schema=schema_correct_external,
        features={"feature_a": "z", "feature_b": "y"},
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"feature_a": "z", "feature_b": "y"}

    # call this again without passing features explicitly (they're already part of the artifact)
    curator = ln.curators.DataFrameCurator(
        artifact,
        schema=schema_correct_external,
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"feature_a": "z", "feature_b": "y"}

    # clean up everything
    inferred_schema = artifact.schemas.all()[0]
    artifact.schemas.remove(inferred_schema.id)
    inferred_schema.delete(permanent=True)
    artifact.delete(permanent=True)
    schema_with_mistake.delete(permanent=True)
    schema_no_external.delete(permanent=True)
    schema_correct_external.delete(permanent=True)
    schema_external.delete(permanent=True)


def test_artifact_save_validates_schema_suffix():
    csv_path = "test_schema_suffix_validation.csv"
    with open(csv_path, "w") as file:
        file.write("a,b\n1,2\n")

    # maximal_set keeps this hash distinct from the shared valid_features schema
    # once suffix is cleared. That schema already has itype=Feature and no suffix.
    schema = ln.Schema(itype=ln.Feature, suffix=".csv", maximal_set=True).save()
    artifact = ln.Artifact(csv_path, key=csv_path, schema=schema).save()

    schema.suffix = ".h5ad"
    schema.save()
    with pytest.raises(ln.errors.ValidationError) as error:
        ln.Artifact(csv_path, key=csv_path, schema=schema).save()
    assert (
        "Artifact not validated. Artifact suffix '.csv' does not match schema "
        "suffix '.h5ad'." in error.exconly()
    )

    schema.suffix = None
    schema.save()
    artifact_none_suffix = ln.Artifact(csv_path, key=csv_path, schema=schema).save()

    for saved_artifact in [artifact, artifact_none_suffix]:
        saved_artifact.delete(permanent=True)

    schema.delete(permanent=True)
    Path(csv_path).unlink(missing_ok=True)


def test_artifact_constructor_does_not_attach_unvalidated_schema(tmp_path: Path):
    gene = ln.Feature(name="attach_schema_gene", dtype=str).save()
    score = ln.Feature(name="attach_schema_score", dtype=float).save()
    schema = ln.Schema(
        name="attach_schema_needs_gene_score",
        features=[gene, score],
    ).save()
    # maximal_set keeps this hash distinct from `schema`, and the missing name
    # exercises the unlabeled error
    unnamed = ln.Schema(features=[gene, score], maximal_set=True).save()
    bad = pd.DataFrame({"attach_schema_gene": ["a", "b"], "other": [1, 2]})
    parquet = tmp_path / "bad.parquet"
    csv = tmp_path / "bad.csv"
    tsv = tmp_path / "bad.tsv"
    text = tmp_path / "notes.txt"
    json_path = tmp_path / "data.json"
    bad.to_parquet(parquet)
    bad.to_csv(csv, index=False)
    bad.to_csv(tsv, sep="\t", index=False)
    text.write_text("not a table")
    json_path.write_text('{"x": 1}')

    with pytest.raises(ln.errors.ValidationError):
        ln.Artifact.from_dataframe(bad, key="attach-schema/a.parquet", schema=schema)

    for path, key in (
        (parquet, "attach-schema/b.parquet"),
        (csv, "attach-schema/c.csv"),
    ):
        with pytest.raises(ln.errors.InvalidArgument) as error:
            ln.Artifact(path, key=key, schema=schema)
        assert "from_dataframe" in error.exconly()
        assert "attach_schema_needs_gene_score" in error.exconly()
        assert not ln.Artifact.filter(key=key).exists()

    # `.tsv` does not get otype="DataFrame", so the tabular-suffix branch runs
    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.Artifact(tsv, key="attach-schema/c.tsv", schema=unnamed)
    assert "from_dataframe" in error.exconly()
    assert "Schema '" not in error.exconly()

    for path, key in (
        (text, "attach-schema/d.txt"),
        (json_path, "attach-schema/e.json"),
    ):
        with pytest.raises(ln.errors.ValidationError) as error:
            ln.Artifact(path, key=key, schema=schema)
        assert "no feature values were provided" in error.exconly()
        assert not ln.Artifact.filter(key=key).exists()

    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.Artifact(text, key="attach-schema/d.txt", otype="FASTQ", schema=schema)
    assert "otype 'FASTQ'" in error.exconly()
    assert "__external__" in error.exconly()

    # same hash as `parquet`, so this is the deduplication path
    saved = ln.Artifact(parquet, key="attach-schema/existing.parquet").save()
    assert saved.schema is None
    with pytest.raises(ln.errors.InvalidArgument):
        ln.Artifact(parquet, key="attach-schema/existing.parquet", schema=schema)
    assert ln.Artifact.get(id=saved.id).schema_id is None
    assert schema.validated_artifacts.count() == 0

    with pytest.raises(ln.errors.ValidationError):
        ln.Artifact(
            text,
            key="attach-schema/d.txt",
            schema=schema,
            features={"attach_schema_gene": "a"},
        )

    annotated = ln.Artifact(
        text,
        key="attach-schema/d.txt",
        schema=schema,
        features={"attach_schema_gene": "a", "attach_schema_score": 1.0},
    ).save()
    assert annotated.schema == schema
    assert schema.validated_artifacts.count() == 1

    # existing artifact supplies the feature values
    again = ln.Artifact(text, key="attach-schema/d.txt", schema=schema).save()
    assert again.uid == annotated.uid
    assert again.features.get_values() == {
        "attach_schema_gene": "a",
        "attach_schema_score": 1.0,
    }

    annotated.delete(permanent=True)
    saved.delete(permanent=True)
    schema.delete(permanent=True)
    unnamed.delete(permanent=True)
    gene.delete(permanent=True)
    score.delete(permanent=True)


def test_external_only_schema_on_any_artifact(
    two_external_features: tuple[ln.Feature, ln.Feature],
    tmp_path: Path,
):
    feat_a, feat_b = two_external_features
    # maximal_set rejects features that are not in the schema
    external = ln.Schema(features=[feat_a, feat_b], maximal_set=True).save()
    schema = ln.Schema(
        name="external_only_any_artifact",
        slots={"__external__": external},
    ).save()
    assert schema.otype is None

    fastq = tmp_path / "reads.fastq.gz"
    fastq.write_bytes(b"@read1\nACGT\n+\n!!!!\n")
    key = "seq/reads.fastq.gz"
    values = {"feature_a": "x", "feature_b": "y"}

    with pytest.raises(ln.errors.ValidationError) as error:
        ln.Artifact(fastq, key=key, otype="FASTQ", schema=schema)
    assert (
        "External features slot is defined in schema but no external features were provided."
        in error.exconly()
    )
    with pytest.raises(ln.errors.ValidationError):
        ln.Artifact(
            fastq,
            key=key,
            otype="FASTQ",
            features={"feature_a": 1, "feature_b": "y"},
            schema=schema,
        )
    extra = ln.Feature(name="not_in_external_schema", dtype=str).save()
    with pytest.raises(ln.errors.ValidationError):
        ln.Artifact(
            fastq,
            key=key,
            otype="FASTQ",
            features={**values, "not_in_external_schema": "z"},
            schema=schema,
        )

    artifact = ln.Artifact(
        fastq,
        key=key,
        otype="FASTQ",
        features=values,
        schema=schema,
    ).save()
    assert artifact.otype == "FASTQ"
    assert artifact.schema == schema
    assert artifact.features.get_values() == values
    described = artifact.features.describe(return_str=True)
    assert "External features" in described
    assert "feature_a" in described

    with pytest.raises(ln.errors.ValidationError):
        artifact.features.set_values({"feature_a": 1, "feature_b": "y"})
    unknown_values: dict[str | ln.Feature, str] = {**values, extra: "z"}
    with pytest.raises(ln.errors.ValidationError) as error:
        artifact.features.set_values(unknown_values)
    assert "not in the provided schema" in error.exconly()
    artifact.features.set_values({"feature_a": "z", "feature_b": "y"})
    assert artifact.features.get_values()["feature_a"] == "z"
    with pytest.raises(ValueError) as error:
        artifact.features.remove_values("feature_a", value="z")
    assert "Cannot remove values if artifact has external schema." in error.exconly()

    again = ln.Artifact(fastq, key=key, otype="FASTQ", schema=schema).save()
    assert again.uid == artifact.uid
    assert again.features.get_values() == {"feature_a": "z", "feature_b": "y"}

    folder = tmp_path / "sample_folder"
    folder.mkdir()
    (folder / "notes.txt").write_text("later")
    folder_artifact = ln.Artifact(
        folder,
        key="seq/sample_folder",
        features=values,
        schema=schema,
    ).save()
    assert folder_artifact.otype is None
    assert folder_artifact.schema == schema
    assert folder_artifact.features.get_values() == values

    content_feature = ln.Feature(name="obs_col_block", dtype=str).save()
    obs_schema = ln.Schema(features=[content_feature]).save()
    content_schema = ln.Schema(
        otype="AnnData",
        slots={"obs": obs_schema, "__external__": external},
    ).save()
    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.Artifact(fastq, key="seq/blocked.fastq.gz", schema=content_schema)
    assert "from_anndata" in error.exconly()

    column_schema = ln.Schema(
        features=[feat_a],
        slots={"__external__": external},
        otype="DataFrame",
    ).save()
    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.Artifact(
            fastq,
            key="seq/blocked-df.fastq.gz",
            features=values,
            schema=column_schema,
        )
    assert "from_dataframe" in error.exconly()

    folder_artifact.delete(permanent=True)
    artifact.delete(permanent=True)
    column_schema.delete(permanent=True)
    content_schema.delete(permanent=True)
    obs_schema.delete(permanent=True)
    content_feature.delete(permanent=True)
    schema.delete(permanent=True)
    external.delete(permanent=True)
    extra.delete(permanent=True)


def test_from_anndata_with_external_schema(
    two_external_features: tuple[ln.Feature, ln.Feature],
):
    import anndata as ad
    import numpy as np

    feat_a, feat_b = two_external_features
    obs_feature = ln.Feature(name="obs_note_external_slot", dtype=str).save()
    obs_schema = ln.Schema(features=[obs_feature]).save()
    external = ln.Schema(features=[feat_a, feat_b]).save()
    schema = ln.Schema(
        name="anndata_with_external_slot",
        otype="AnnData",
        slots={"obs": obs_schema, "__external__": external},
    ).save()
    adata = ad.AnnData(
        np.array([[1.0, 2.0], [3.0, 4.0]]),
        obs=pd.DataFrame({"obs_note_external_slot": ["a", "b"]}),
    )

    with pytest.raises(ln.errors.ValidationError) as error:
        ln.Artifact.from_anndata(adata, key="external_slot.h5ad", schema=schema)
    assert (
        "External features slot is defined in schema but no external features were provided."
        in error.exconly()
    )

    artifact = ln.Artifact.from_anndata(
        adata,
        key="external_slot.h5ad",
        schema=schema,
        features={"feature_a": "x", "feature_b": "y"},
    ).save()
    assert artifact.features.get_values()["feature_a"] == "x"
    assert artifact.features.get_values()["feature_b"] == "y"
    described = artifact.features.describe(return_str=True)
    assert "External features" in described
    assert "obs" in described

    inferred = list(artifact.schemas.all())
    for inferred_schema in inferred:
        artifact.schemas.remove(inferred_schema)
        if inferred_schema.id not in {schema.id, obs_schema.id, external.id}:
            inferred_schema.delete(permanent=True)
    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    obs_schema.delete(permanent=True)
    external.delete(permanent=True)
    obs_feature.delete(permanent=True)
