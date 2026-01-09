import os

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from django.db.utils import IntegrityError
from lamindb.errors import FieldValidationError, InvalidArgument, ValidationError
from lamindb.models.schema import (
    get_related_name,
    migrate_auxiliary_fields_postgres,
    validate_features,
)


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
    bt.Gene.filter(symbol__in=gene_symbols).delete(permanent=True)
    with pytest.raises(ValidationError) as error:
        schema = ln.Schema.from_values(
            gene_symbols, bt.Gene.symbol, dtype=int, organism="human"
        )
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: These values could not be validated:"
    )
    ln.save(bt.Gene.from_values(gene_symbols, "symbol", organism="human"))
    schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol, organism="human")
    # below should be a queryset and not a list
    assert set(schema.members) == set(
        bt.Gene.from_values(gene_symbols, "symbol", organism="human")
    )
    assert schema.dtype == "num"  # this is NUMBER_TYPE
    schema = ln.Schema.from_values(
        gene_symbols, bt.Gene.symbol, dtype=int, organism="human"
    )
    assert schema._state.adding
    assert schema.dtype == "int"
    assert schema.itype == "bionty.Gene"
    schema.save()
    assert set(schema.members) == set(schema.genes.all())
    id = schema.id
    # test that the schema is retrieved from the database
    # in case it already exists
    schema = ln.Schema.from_values(
        gene_symbols, bt.Gene.symbol, dtype=int, organism="human"
    )
    assert not schema._state.adding
    assert id == schema.id
    schema.delete(permanent=True)

    # edge cases
    with pytest.raises(ValueError):
        schema = ln.Schema.from_values([])
    with pytest.raises(TypeError):
        ln.Schema.from_values(["a"], field="name")
    with pytest.raises(ValidationError):
        schema = ln.Schema.from_values(
            ["weird_name"], field=ln.Feature.name, dtype="float"
        )


def test_schema_from_records(df):
    features = ln.Feature.from_dataframe(df)
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
    schema.delete(permanent=True)

    # edge case
    with pytest.raises(ValueError):
        positional_arg = 1
        ln.Schema(features, positional_arg)


def test_schema_from_df(df):
    # test using type
    human = bt.Organism.from_source(name="human").save()
    genes = [bt.Gene(symbol=name, organism=human) for name in df.columns]
    ln.save(genes)
    with pytest.raises(ValueError) as error:
        ln.Schema.from_dataframe(df, field=bt.Gene.symbol)
    assert error.exconly().startswith("ValueError: data types are heterogeneous:")
    schema = ln.Schema.from_dataframe(df[["feat1", "feat2"]], field=bt.Gene.symbol)
    for gene in genes:
        gene.delete(permanent=True)

    # now for the features registry
    features = ln.Feature.from_dataframe(df)
    ln.save(features)
    schema = ln.Schema.from_dataframe(df).save()
    assert schema.dtype is None
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


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
    transform = ln.Transform(key="test").save()
    # This is just a type check
    with pytest.raises(TypeError) as error:
        validate_features([transform, ln.Run(transform)])
    assert error.exconly() == "TypeError: schema can only contain a single type"
    transform.delete(permanent=True)


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
    feature.delete(permanent=True)


@pytest.fixture(scope="module")
def mini_immuno_schema_flexible():
    schema = ln.examples.datasets.mini_immuno.define_mini_immuno_schema_flexible()

    yield schema

    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    bt.Gene.filter().delete(permanent=True)
    ln.Record.filter(type__isnull=False).delete(permanent=True)
    ln.Record.filter().delete(permanent=True)
    bt.CellType.filter().delete(permanent=True)


def test_schema_update_implicit_through_name_equality(
    mini_immuno_schema_flexible: ln.Schema,
    ccaplog,
):
    df = pd.DataFrame({"a": [1]})
    artifact = ln.Artifact.from_dataframe(df, key="test_artifact.parquet").save()
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
    ln.examples.datasets.mini_immuno.define_mini_immuno_schema_flexible()

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
    artifact = ln.Artifact.from_dataframe(df, key="test_artifact.parquet").save()
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
    feature_to_add.delete(permanent=True)

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

    # change coerce (formerly auxiliary field, now Django field) --------------------------------

    assert not mini_immuno_schema_flexible.coerce
    mini_immuno_schema_flexible.coerce = True
    mini_immuno_schema_flexible.save()
    assert mini_immuno_schema_flexible.hash != orig_hash
    assert ccaplog.text.count(warning_message) == 5

    # restore original setting
    mini_immuno_schema_flexible.coerce = False
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
    index_feature.delete(permanent=True)

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


def test_schema_mutations_feature_removal(
    mini_immuno_schema_flexible: ln.Schema, ccaplog
):
    feature1 = ln.Feature.get(name="perturbation")
    feature2 = ln.Feature.get(name="cell_type_by_model")
    dummy_artifact = ln.Artifact(".gitignore", key=".gitignore").save()
    # define the schema the first time
    schema = ln.Schema(name="My test schema X", features=[feature1, feature2]).save()
    assert schema.features.count() == 2
    dummy_artifact.schema = schema  # pretend artifact was validated with this schema
    dummy_artifact.save()
    # define the schema the first time
    schema1 = ln.Schema(name="My test schema X", features=[feature2]).save()
    # retrieves same schema because of name equality
    assert ccaplog.text.count("you're removing these features:") == 1
    assert (
        ccaplog.text.count("you updated the schema hash and might invalidate datasets")
        == 1
    )
    assert schema1 == schema
    assert schema1.features.count() == 1
    dummy_artifact.delete(permanent=True)
    schema.delete(permanent=True)


def test_schema_add_remove_optional_features(mini_immuno_schema_flexible: ln.Schema):
    schema = mini_immuno_schema_flexible
    initial_hash = schema.hash
    feature_project = ln.Feature(name="project", dtype=ln.Project).save()
    schema.add_optional_features([feature_project])
    assert schema.hash != initial_hash
    schema.remove_optional_features([feature_project])
    assert schema.hash == initial_hash


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

    with pytest.raises(InvalidArgument) as error:
        ln.Schema(
            name="mini_immuno_anndata_schema",
            slots={"obs": obs_schema, "var": var_schema},
        ).save()
    assert str(error.value) == "Please pass otype != None for composite schemas"

    anndata_schema = ln.Schema(
        name="mini_immuno_anndata_schema",
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
    with pytest.raises(IntegrityError) as error:
        anndata_schema.components.add(var_schema2, through_defaults={"slot": "var"})
    assert "unique" in str(error.value).lower()

    anndata_schema.delete(permanent=True)
    var_schema2.delete(permanent=True)
    var_schema.delete(permanent=True)


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


def test_schema_recovery_based_on_hash(mini_immuno_schema_flexible: ln.Schema):
    feature1 = ln.Feature.get(name="perturbation")
    feature2 = ln.Feature.get(name="cell_type_by_model")
    schema = ln.Schema(features=[feature1, feature2]).save()
    schema2 = ln.Schema(features=[feature1, feature2])
    assert schema == schema2
    schema.delete()
    schema2 = ln.Schema(features=[feature1, feature2])
    assert schema != schema2
    schema.delete(permanent=True)


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
        "a=num",
        "c=True",
        "d=False",
        "e=False",
        "l=GPZ-TzvKRhdC1PQAhlFiow",
    ]
    assert schema.name == "anndata_ensembl_gene_ids_and_valid_features_in_obs"
    assert schema.itype is None
    assert schema.hash == "aqGWHvyY49W_PHELUMiBMw"

    # test the convenience function
    schema = ln.examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs()
    assert schema.uid == "0000000000000002"
    assert schema.name == "anndata_ensembl_gene_ids_and_valid_features_in_obs"
    assert schema.itype is None
    assert schema.hash == "aqGWHvyY49W_PHELUMiBMw"
    varT_schema = schema.slots["var.T"]
    assert varT_schema.uid == "0000000000000001"
    assert varT_schema.name == "valid_ensembl_gene_ids"
    assert varT_schema.itype == "bionty.Gene.ensembl_gene_id"
    assert varT_schema.hash == "1gocc_TJ1RU2bMwDRK-WUA"

    schema.delete(permanent=True)


def test_schema_already_saved_aux():
    """When attempting to save a Schema that was already saved before which populated `_aux` fields,
    we expect the Schema to be returned with the same `_aux` fields.

    Test for https://github.com/laminlabs/lamindb/issues/2887
    """
    var_schema = ln.Schema(
        name="test var",
        index=ln.Feature(
            name="var_index",
            dtype=bt.Gene.ensembl_gene_id,
            cat_filters={
                "source": bt.Source.get(
                    entity="bionty.Gene", currently_used=True, organism="human"
                )
            },
        ).save(),
        itype=ln.Feature,
        dtype="DataFrame",
        minimal_set=True,
        coerce=True,
    ).save()

    schema = ln.Schema(
        name="AnnData schema",
        otype="AnnData",
        minimal_set=True,
        coerce=True,
        slots={"var": var_schema},
    ).save()

    # _aux["af"] now only contains key "3" (index_feature_uid) since coerce and flexible are Django fields
    assert len(schema.slots["var"]._aux["af"].keys()) == 1
    assert "3" in schema.slots["var"]._aux["af"]  # index_feature_uid
    # coerce and flexible are now proper Django fields
    assert schema.slots["var"].coerce is True
    assert schema.slots["var"].flexible is False

    # Attempting to save the same schema again should return the Schema with the same fields
    var_schema_2 = ln.Schema(
        name="test var",
        index=ln.Feature(
            name="var_index",
            dtype=bt.Gene.ensembl_gene_id,
            cat_filters={
                "source": bt.Source.get(
                    entity="bionty.Gene", currently_used=True, organism="human"
                )
            },
        ).save(),
        itype=ln.Feature,
        dtype="DataFrame",
        minimal_set=True,
        coerce=True,
    ).save()

    schema_2 = ln.Schema(
        name="AnnData schema",
        otype="AnnData",
        minimal_set=True,
        coerce=True,
        slots={"var": var_schema_2},
    ).save()

    assert len(schema.slots["var"]._aux["af"].keys()) == 1
    assert schema.slots["var"]._aux == schema_2.slots["var"]._aux
    assert schema.slots["var"].coerce == schema_2.slots["var"].coerce
    assert schema.slots["var"].flexible == schema_2.slots["var"].flexible

    schema_2.delete(permanent=True)
    schema.delete(permanent=True)


def test_schema_not_saved_describe():
    schema = ln.Schema(name="NotSavedSchema", is_type=True)
    with pytest.raises(ValueError) as e:
        schema.describe()
    assert "Schema must be saved before describing" in str(e.value)


def test_schema_is_type():
    Sample = ln.Schema(name="Sample", is_type=True).save()
    assert Sample.hash is None
    BioSample = ln.Schema(name="BioSample", is_type=True, type=Sample).save()
    assert BioSample.hash is None
    assert BioSample.type == Sample
    assert BioSample.is_type
    # create a schema without any features or slots or itype or is_type=True
    with pytest.raises(InvalidArgument) as e:
        ln.Schema(name="TechSample", type=Sample)
    assert "Please pass features or slots or itype or set is_type=True" in str(e.value)
    # clean up
    BioSample.delete(permanent=True)
    Sample.delete(permanent=True)


# see test_component_composite in test_transform.py
def test_composite_component():
    composite = ln.Schema(name="composite", itype=ln.Feature).save()
    component1 = ln.Schema(name="component1", itype=bt.CellType).save()
    component2 = ln.Schema(name="component2", itype=bt.CellMarker).save()
    composite.components.add(component1, through_defaults={"slot": "slot1"})
    composite.components.add(component2, through_defaults={"slot": "slot2"})

    assert len(composite.components.all()) == 2
    assert composite.links_component.count() == 2
    assert set(composite.links_component.all().to_list("slot")) == {"slot1", "slot2"}
    assert composite.links_component.first().composite == composite
    assert composite.composites.count() == 0
    assert composite.links_composite.count() == 0

    ln.models.SchemaComponent.filter(composite=composite).delete(permanent=True)

    link = ln.models.SchemaComponent(
        composite=composite, component=component1, slot="var"
    ).save()
    assert link in composite.links_component.all()
    assert link in component1.links_composite.all()
    assert link.slot == "var"

    composite.delete(permanent=True)
    component1.delete(permanent=True)
    component2.delete(permanent=True)

    assert ln.models.SchemaComponent.filter().count() == 0


@pytest.mark.skipif(
    os.getenv("LAMINDB_TEST_DB_VENDOR") != "postgresql",
    reason="PostgreSQL-specific migration test",
)
def test_migrate_auxiliary_fields_postgres():
    """Test PostgreSQL migration of auxiliary fields for all models.

    This test verifies that migrate_auxiliary_fields_postgres correctly migrates:

    **Artifact:**
    - _storage_completed from _aux['af']['0']

    **Run:**
    - cli_args from _aux['af']['0']

    **Feature:**
    - default_value from _aux['af']['0']
    - nullable from _aux['af']['1'] (default: True)
    - coerce from _aux['af']['2'] (default: False)
    - For type features, all values are set to NULL

    **Schema:**
    - coerce from _aux['af']['0']
    - flexible from _aux['af']['2'] (or computes from n_members)
    - Converts negative n_members to NULL
    - For type schemas, all values are set to NULL
    - Preserves '1' (optionals) and '3' (index_feature_uid) in _aux
    """
    from django.db import connection

    # === Setup test data ===

    # Create a Transform and Run for testing
    transform = ln.Transform(key="test_migration_transform").save()
    run = ln.Run(transform=transform).save()

    # Create an Artifact for testing
    artifact = ln.Artifact(".gitignore", key="test_migration_artifact").save()

    # Create Features for testing (type and regular)
    type_feature = ln.Feature(
        name="TestMigrationTypeFeat", dtype=str, is_type=True
    ).save()
    regular_feature = ln.Feature(name="test_migration_regular_feat", dtype=str).save()

    # Create Schemas for testing (type and regular)
    type_schema = ln.Schema(name="TestMigrationTypeSchema", is_type=True).save()
    feature_for_schema1 = ln.Feature(
        name="test_migration_schema_feat1", dtype=str
    ).save()
    feature_for_schema2 = ln.Feature(
        name="test_migration_schema_feat2", dtype=str
    ).save()
    regular_schema = ln.Schema(
        name="TestMigrationRegularSchema",
        features=[feature_for_schema1, feature_for_schema2],
        coerce=True,
        flexible=True,
    ).save()

    # === Set old-style _aux data to simulate pre-migration state ===
    with connection.cursor() as cursor:
        # Artifact: set _aux with af containing _storage_completed value
        cursor.execute(
            """
            UPDATE lamindb_artifact
            SET _aux = '{"af": {"0": true}}'::jsonb,
                _storage_completed = NULL
            WHERE id = %s
            """,
            [artifact.id],
        )

        # Run: set _aux with af containing cli_args value
        cursor.execute(
            """
            UPDATE lamindb_run
            SET _aux = '{"af": {"0": "--verbose --debug"}}'::jsonb,
                cli_args = NULL
            WHERE id = %s
            """,
            [run.id],
        )

        # Feature (type): set _aux with af keys that should result in NULL values
        cursor.execute(
            """
            UPDATE lamindb_feature
            SET _aux = '{"af": {"0": "default_val", "1": false, "2": true}}'::jsonb,
                default_value = NULL,
                nullable = NULL,
                coerce = NULL
            WHERE id = %s
            """,
            [type_feature.id],
        )

        # Feature (regular): set _aux with af keys for migration
        cursor.execute(
            """
            UPDATE lamindb_feature
            SET _aux = '{"af": {"0": "my_default", "1": false, "2": true}}'::jsonb,
                default_value = NULL,
                nullable = NULL,
                coerce = NULL
            WHERE id = %s
            """,
            [regular_feature.id],
        )

        # Schema (type): set _aux with af keys that should be cleaned
        cursor.execute(
            """
            UPDATE lamindb_schema
            SET _aux = '{"af": {"0": true, "2": false}}'::jsonb,
                coerce = NULL,
                flexible = NULL
            WHERE id = %s
            """,
            [type_schema.id],
        )

        # Schema (regular): set _aux with af keys including optionals (key "1")
        cursor.execute(
            """
            UPDATE lamindb_schema
            SET _aux = '{"af": {"0": true, "1": ["uid1", "uid2"], "2": true}}'::jsonb,
                coerce = NULL,
                flexible = NULL
            WHERE id = %s
            """,
            [regular_schema.id],
        )

    # === Run the migration function ===
    with connection.schema_editor() as schema_editor:
        migrate_auxiliary_fields_postgres(schema_editor)

    # === Refresh all objects from database ===
    artifact.refresh_from_db()
    run.refresh_from_db()
    type_feature.refresh_from_db()
    regular_feature.refresh_from_db()
    type_schema.refresh_from_db()
    regular_schema.refresh_from_db()

    # === Verify Artifact migration ===
    assert artifact._storage_completed is True  # from _aux['af']['0']
    # _aux should have 'af' removed (was only key)
    assert artifact._aux is None or "af" not in artifact._aux

    # === Verify Run migration ===
    assert run.cli_args == "--verbose --debug"  # from _aux['af']['0']
    # _aux should have 'af' removed
    assert run._aux is None or "af" not in run._aux

    # === Verify Feature (type) migration ===
    # Type features should have all values set to NULL
    assert type_feature.default_value is None
    assert type_feature.nullable is None
    assert type_feature.coerce is None
    # _aux should have 'af' removed
    assert type_feature._aux is None or "af" not in type_feature._aux

    # === Verify Feature (regular) migration ===
    assert regular_feature.default_value == "my_default"  # from _aux['af']['0']
    assert regular_feature.nullable is False  # from _aux['af']['1']
    assert regular_feature.coerce is True  # from _aux['af']['2']
    # _aux should have 'af' removed
    assert regular_feature._aux is None or "af" not in regular_feature._aux

    # === Verify Schema (type) migration ===
    assert type_schema.coerce is None
    assert type_schema.flexible is None
    assert type_schema.n_members is None
    # _aux should either be None or not have '0' and '2' keys in 'af'
    if type_schema._aux is not None and "af" in type_schema._aux:
        assert "0" not in type_schema._aux["af"]
        assert "2" not in type_schema._aux["af"]

    # === Verify Schema (regular) migration ===
    assert regular_schema.coerce is True  # from _aux['af']['0']
    assert regular_schema.flexible is True  # from _aux['af']['2']
    # _aux should preserve key '1' (optionals)
    assert regular_schema._aux is not None
    assert "af" in regular_schema._aux
    assert "1" in regular_schema._aux["af"]
    assert regular_schema._aux["af"]["1"] == ["uid1", "uid2"]
    # Keys '0' and '2' should be removed
    assert "0" not in regular_schema._aux["af"]
    assert "2" not in regular_schema._aux["af"]

    # === Clean up ===
    regular_schema.delete(permanent=True)
    type_schema.delete(permanent=True)
    feature_for_schema1.delete(permanent=True)
    feature_for_schema2.delete(permanent=True)
    regular_feature.delete(permanent=True)
    type_feature.delete(permanent=True)
    artifact.delete(permanent=True)
    run.delete(permanent=True)
    transform.delete(permanent=True)
