import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from django.db.utils import IntegrityError
from lamindb.errors import ValidationError
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
    with pytest.raises(ValueError):
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
def small_dataset1_schema():
    # define labels
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
    ln.ULabel(name="DMSO", type=perturbation).save()
    ln.ULabel(name="IFNG", type=perturbation).save()
    bt.CellType.from_source(name="B cell").save()
    bt.CellType.from_source(name="T cell").save()

    # in next iteration for attrs
    # ln.Feature(name="temperature", dtype="float").save()
    # ln.Feature(name="study", dtype="cat[ULabel]").save()
    # ln.Feature(name="date_of_study", dtype="date").save()
    # ln.Feature(name="study_note", dtype="str").save()

    # define schema
    schema = ln.Schema(
        name="small_dataset1_obs_level_metadata",
        features=[
            ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
            ln.Feature(name="sample_note", dtype=str).save(),
            ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
            ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
        ],
    ).save()

    yield schema

    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter(type__isnull=False).delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()


def test_schema_recreation_with_same_name_different_hash(
    small_dataset1_schema: ln.Schema,
):
    try:
        ln.Schema(
            name="small_dataset1_obs_level_metadata",
            features=[
                ln.Feature.get(name="perturbation"),
                ln.Feature.get(name="sample_note"),
            ],
        ).save()
    except ValueError as error:
        assert str(error).startswith("Schema name is already in use by schema with uid")


def test_schema_components(small_dataset1_schema: ln.Schema):
    obs_schema = small_dataset1_schema
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
            components={"obs": obs_schema, "var": var_schema},
        ).save()
    except ln.errors.InvalidArgument:
        assert (
            str(ln.errors.InvalidArgument)
            == "Please pass otype != None for composite schemas"
        )

    anndata_schema = ln.Schema(
        name="small_dataset1_anndata_schema",
        otype="AnnData",
        components={"obs": obs_schema, "var": var_schema},
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
