import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from django.db.utils import IntegrityError
from lamindb.errors import FieldValidationError, ValidationError
from lamindb.models.schema import get_related_name, validate_features


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame(
        {
            "feat1": [1, 2, 3],
            "feat2": [3, 4, 5],
            "feat3": ["cond1", "cond2", "cond2"],
            "feat4": ["id1", "id2", "id3"],
        }
    )


def test_schema_from_values():
    gene_symbols = ["TCF7", "MYC"]
    bt.settings.organism = "human"
    bt.Gene.filter(symbol__in=gene_symbols).all().delete()
    with pytest.raises(ValidationError) as error:
        schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: These values could not be validated:"
    )
    ln.save(bt.Gene.from_values(gene_symbols, "symbol"))
    schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol)
    # below should be a queryset and not a list
    assert set(schema.members) == set(bt.Gene.from_values(gene_symbols, "symbol"))
    assert schema.dtype == "num"  # this is NUMBER_TYPE
    schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert schema._state.adding
    assert schema.dtype == "int"
    assert schema.itype == "bionty.Gene"
    schema.save()
    assert set(schema.members) == set(schema.genes.all())
    id = schema.id
    # test that the schema is retrieved from the database
    # in case it already exists
    schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert not schema._state.adding
    assert id == schema.id
    schema.delete()

    # edge cases
    with pytest.raises(ValueError):
        schema = ln.Schema.from_values([])

    with pytest.raises(TypeError):
        ln.Schema.from_values(["a"], field="name")
    with pytest.raises(ValidationError):
        schema = ln.Schema.from_values(
            ["weird_name"], field=ln.Feature.name, type="float"
        )
    with pytest.raises(ValidationError):
        ln.Schema.from_values([1], field=ln.Feature.name, type="float")

    # return none if no validated features
    with pytest.raises(ValidationError):
        ln.Schema.from_values(["name"], field=ln.Feature.name, type="float")


def test_schema_from_records(df):
    features = ln.Feature.from_df(df)
    with pytest.raises(ValueError) as error:
        schema = ln.Schema(features)
    assert (
        error.exconly()
        == "ValueError: Can only construct feature sets from validated features"
    )

    ln.save(features)
    schema = ln.Schema(features)
    assert schema.id is None
    assert schema._state.adding
    assert schema.dtype is None
    assert schema.itype == "Feature"
    schema.save()
    # test that the schema is retrieved from the database
    # in case it already exists
    schema = ln.Schema(features)
    assert not schema._state.adding
    assert schema.id is not None
    schema.delete()

    # edge case
    with pytest.raises(ValueError):
        positional_arg = 1
        ln.Schema(features, positional_arg)


def test_schema_from_df(df):
    # test using type
    bt.settings.organism = "human"
    genes = [bt.Gene(symbol=name) for name in df.columns]
    ln.save(genes)
    with pytest.raises(ValueError) as error:
        ln.Schema.from_df(df, field=bt.Gene.symbol)
    assert error.exconly().startswith("ValueError: data types are heterogeneous:")
    schema = ln.Schema.from_df(df[["feat1", "feat2"]], field=bt.Gene.symbol)
    for gene in genes:
        gene.delete()

    # now for the features registry
    features = ln.Feature.from_df(df)
    ln.save(features)
    schema = ln.Schema.from_df(df).save()
    assert schema.dtype is None
    ln.Schema.filter().all().delete()
    ln.Feature.filter().all().delete()


def test_get_related_name():
    with pytest.raises(ValueError):
        get_related_name(ln.Transform)


def test_validate_features():
    with pytest.raises(ValueError):
        validate_features([])
    with pytest.raises(TypeError):
        validate_features(["feature"])
    with pytest.raises(TypeError):
        validate_features({"feature"})
    transform = ln.Transform(key="test")
    transform.save()
    # This is just a type check
    with pytest.raises(TypeError) as error:
        validate_features([transform, ln.Run(transform)])
    assert error.exconly() == "TypeError: schema can only contain a single type"
    transform.delete()


def test_kwargs():
    with pytest.raises(FieldValidationError):
        ln.Schema(x="1", features=[])


def test_edge_cases():
    feature = ln.Feature(name="rna", dtype="float")
    ln.save([feature])
    with pytest.raises(ValueError) as error:
        ln.Schema(feature)
    assert (
        error.exconly()
        == "ValueError: Please pass a ListLike of features, not a single feature"
    )
    feature.delete()


@pytest.fixture(scope="module")
def mini_immuno_schema_flexible():
    schema = ln.core.datasets.mini_immuno.define_mini_immuno_schema_flexible()

    yield schema

    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter(type__isnull=False).delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()


def test_schema_update_implicit_through_name_equality(
    mini_immuno_schema_flexible: ln.Schema,
    ccaplog,
):
    df = pd.DataFrame({"a": [1]})
    artifact = ln.Artifact.from_df(df, key="test_artifact.parquet").save()
    artifact.schema = mini_immuno_schema_flexible
    artifact.save()

    orig_hash = mini_immuno_schema_flexible.hash
    warning_message = "you updated the schema hash and might invalidate datasets that were previously validated with this schema:"

    # different numbers of features -------------------------------------------

    schema = ln.Schema(
        name="Mini immuno schema",
        features=[
            ln.Feature.get(name="perturbation"),
            ln.Feature.get(name="donor"),
        ],
    ).save()

    assert schema.hash != orig_hash
    assert ccaplog.text.count(warning_message) == 1

    # change is flexible (an auxiliary field) --------------------------------

    schema = ln.Schema(
        name="Mini immuno schema",
        features=[
            ln.Feature.get(name="perturbation"),
            ln.Feature.get(name="cell_type_by_model"),
            ln.Feature.get(name="assay_oid"),
            ln.Feature.get(name="donor"),
            ln.Feature.get(name="concentration"),
            ln.Feature.get(name="treatment_time_h"),
        ],
        flexible=True,
    ).save()

    assert schema.hash == orig_hash  # restored original hash
    assert ccaplog.text.count(warning_message) == 2  # warning raised

    schema = ln.Schema(
        name="Mini immuno schema",
        features=[
            ln.Feature.get(name="perturbation"),
            ln.Feature.get(name="cell_type_by_model"),
            ln.Feature.get(name="assay_oid"),
            ln.Feature.get(name="donor"),
            ln.Feature.get(name="concentration"),
            ln.Feature.get(name="treatment_time_h"),
        ],
        flexible=False,
    ).save()

    assert schema.hash != orig_hash
    assert ccaplog.text.count(warning_message) == 3  # warning raised
    ln.core.datasets.mini_immuno.define_mini_immuno_schema_flexible()

    artifact.delete(permanent=True)

    # restore original hash  --------------------------------

    schema = ln.Schema(
        name="Mini immuno schema",
        features=[
            ln.Feature.get(name="perturbation"),
            ln.Feature.get(name="cell_type_by_model"),
            ln.Feature.get(name="assay_oid"),
            ln.Feature.get(name="donor"),
            ln.Feature.get(name="concentration"),
            ln.Feature.get(name="treatment_time_h"),
        ],
        flexible=True,
    ).save()

    assert schema.hash == orig_hash  # restored original hash


def test_schema_update(
    mini_immuno_schema_flexible: ln.Schema,
    ccaplog,
):
    df = pd.DataFrame({"a": [1]})
    artifact = ln.Artifact.from_df(df, key="test_artifact.parquet").save()
    artifact.schema = mini_immuno_schema_flexible
    artifact.save()

    # store original hash

    orig_hash = mini_immuno_schema_flexible.hash
    warning_message = "you updated the schema hash and might invalidate datasets that were previously validated with this schema:"

    # add a feature -------------------------------------------

    feature_to_add = ln.Feature(name="sample_note", dtype=str).save()
    assert mini_immuno_schema_flexible.n == 6
    mini_immuno_schema_flexible.features.add(feature_to_add)
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash != orig_hash
    assert mini_immuno_schema_flexible.n == 7
    assert ccaplog.text.count(warning_message) == 1

    # remove the feature again
    mini_immuno_schema_flexible.features.remove(feature_to_add)
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash == orig_hash
    assert ccaplog.text.count(warning_message) == 2
    assert mini_immuno_schema_flexible.n == 6
    feature_to_add.delete()

    # change is flexible (an auxiliary field) --------------------------------

    assert mini_immuno_schema_flexible.flexible
    mini_immuno_schema_flexible.flexible = False
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash != orig_hash
    assert ccaplog.text.count(warning_message) == 3

    # restore original setting
    mini_immuno_schema_flexible.flexible = True
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash == orig_hash
    assert ccaplog.text.count(warning_message) == 4

    # change coerce_dtype (an auxiliary field) --------------------------------

    assert not mini_immuno_schema_flexible.coerce_dtype
    mini_immuno_schema_flexible.coerce_dtype = True
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash != orig_hash
    assert ccaplog.text.count(warning_message) == 5

    # restore original setting
    mini_immuno_schema_flexible.coerce_dtype = False
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash == orig_hash
    assert ccaplog.text.count(warning_message) == 6

    # add an index --------------------------------

    index_feature = ln.Feature(name="immuno_sample", dtype=str).save()
    mini_immuno_schema_flexible.index = index_feature
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash != orig_hash
    assert mini_immuno_schema_flexible.n == 7
    assert ccaplog.text.count(warning_message) == 7

    # remove the index
    mini_immuno_schema_flexible.index = None
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.n == 6
    assert mini_immuno_schema_flexible.hash == orig_hash
    assert ccaplog.text.count(warning_message) == 8
    index_feature.delete()

    # make a feature optional --------------------------------

    required_feature = mini_immuno_schema_flexible.features.first()
    mini_immuno_schema_flexible.optionals.add(required_feature)
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash != orig_hash
    assert ccaplog.text.count(warning_message) == 9

    # make it required again
    mini_immuno_schema_flexible.optionals.remove(required_feature)
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash == orig_hash
    assert ccaplog.text.count(warning_message) == 10

    artifact.delete(permanent=True)


def test_schema_components(mini_immuno_schema_flexible: ln.Schema):
    obs_schema = mini_immuno_schema_flexible
    var_schema = ln.Schema(
        name="scRNA_seq_var_schema",
        itype=bt.Gene.ensembl_gene_id,
        dtype="num",
    ).save()

    # test recreation of schema based on name lookup
    var_schema2 = ln.Schema(
        name="scRNA_seq_var_schema",
        itype=bt.Gene.ensembl_gene_id,
        dtype="num",
    ).save()
    assert var_schema == var_schema2

    try:
        ln.Schema(
            name="small_dataset1_anndata_schema",
            otype="AnnData",
            slots={"obs": obs_schema, "var": var_schema},
        ).save()
    except ln.errors.InvalidArgument:
        assert (
            str(ln.errors.InvalidArgument)
            == "Please pass otype != None for composite schemas"
        )

    anndata_schema = ln.Schema(
        name="small_dataset1_anndata_schema",
        otype="AnnData",
        slots={"obs": obs_schema, "var": var_schema},
    ).save()

    var_schema2 = ln.Schema(
        name="symbol_var_schema",
        itype=bt.Gene.symbol,
        dtype="num",
    ).save()
    # try adding another schema under slot "var"
    # we want to trigger the unique constraint on slot
    try:
        anndata_schema.components.add(var_schema2, through_defaults={"slot": "var"})
    except IntegrityError as error:
        assert str(error).startswith("duplicate key value violates unique constraint")

    anndata_schema.delete()
    var_schema2.delete()
    var_schema.delete()


def test_mini_immuno_schema_flexible(mini_immuno_schema_flexible):
    schema = ln.Schema(
        name="Mini immuno schema",
        features=[
            ln.Feature.get(name="perturbation"),
            ln.Feature.get(name="cell_type_by_model"),
            ln.Feature.get(name="assay_oid"),
            ln.Feature.get(name="donor"),
            ln.Feature.get(name="concentration"),
            ln.Feature.get(name="treatment_time_h"),
        ],
        flexible=True,  # _additional_ columns in a dataframe are validated & annotated
    )
    assert schema.name == "Mini immuno schema"
    assert schema.itype == "Feature"
    assert (
        schema._list_for_hashing[:6]
        == [
            "b=Feature",
            "c=True",
            "d=False",
            "e=False",
            "f=True",
            "h=6",
            "j=HASH_OF_FEATURE_UIDS",  # this last hash is not deterministic in a unit test
        ][:6]
    )


def test_schemas_dataframe():
    # test on the Python level after record creation -- no saving!
    schema = ln.Schema(name="valid_features", itype=ln.Feature)
    assert schema.name == "valid_features"
    assert schema.itype == "Feature"
    assert schema._list_for_hashing == [
        "b=Feature",
        "c=True",
        "d=False",
        "e=False",
    ]
    assert schema.hash == "kMi7B_N88uu-YnbTLDU-DA"

    # test the convenience function
    schema = ln.examples.schemas.valid_features()
    assert schema.uid == "0000000000000000"
    assert schema.name == "valid_features"
    assert schema.itype == "Feature"
    assert schema.hash == "kMi7B_N88uu-YnbTLDU-DA"


def test_schemas_anndata():
    # test on the Python level after record creation -- no saving!
    obs_schema = ln.examples.schemas.valid_features()
    varT_schema = ln.Schema(
        name="valid_ensembl_gene_ids", itype=bt.Gene.ensembl_gene_id
    )
    assert varT_schema._list_for_hashing == [
        "a=num",
        "b=bionty.Gene.ensembl_gene_id",
        "c=True",
        "d=False",
        "e=False",
    ]
    assert varT_schema.name == "valid_ensembl_gene_ids"
    assert varT_schema.itype == "bionty.Gene.ensembl_gene_id"
    assert varT_schema.hash == "1gocc_TJ1RU2bMwDRK-WUA"
    schema = ln.Schema(
        name="anndata_ensembl_gene_ids_and_valid_features_in_obs",
        otype="AnnData",
        slots={"obs": obs_schema, "var.T": varT_schema.save()},
    )
    assert schema._list_for_hashing == [
        "1gocc_TJ1RU2bMwDRK-WUA",
        "kMi7B_N88uu-YnbTLDU-DA",
    ]
    assert schema.name == "anndata_ensembl_gene_ids_and_valid_features_in_obs"
    assert schema.itype == "Composite"
    assert schema.hash == "GTxxM36n9tocphLfdbNt9g"

    # test the convenience function
    schema = ln.examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs()
    assert schema.uid == "0000000000000002"
    assert schema.name == "anndata_ensembl_gene_ids_and_valid_features_in_obs"
    assert schema.itype == "Composite"
    assert schema.hash == "GTxxM36n9tocphLfdbNt9g"
    varT_schema = schema.slots["var.T"]
    assert varT_schema.uid == "0000000000000001"
    assert varT_schema.name == "valid_ensembl_gene_ids"
    assert varT_schema.itype == "bionty.Gene.ensembl_gene_id"
    assert varT_schema.hash == "1gocc_TJ1RU2bMwDRK-WUA"
    schema.delete()
