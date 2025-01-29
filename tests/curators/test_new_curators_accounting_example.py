"""Test suite for accounting on bank transactions."""

import datetime

import lamindb as ln
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def currency_labels():
    # Create currency type and labels
    currency_type = ln.ULabel(name="Currency", is_type=True).save()
    usd = ln.ULabel(name="USD", type=currency_type).save()
    eur = ln.ULabel(name="EUR", type=currency_type).save()

    yield currency_type, usd, eur

    # Cleanup
    eur.delete()
    usd.delete()
    currency_type.delete()


@pytest.fixture(scope="module")
def transactions_features():
    # Create features
    currency_feature = ln.Feature(
        name="currency_name", dtype="cat[ULabel[Currency].name]"
    ).save()
    date_feature = ln.Feature(name="date", dtype="date").save()

    transaction_type = ln.Feature(name="Transaction", is_type=True).save()
    amount_usd = ln.Feature(
        name="transaction_amount_usd_cent", dtype=int, type=transaction_type
    ).save()
    amount_eur = ln.Feature(
        name="transaction_amount_eur_cent", dtype=int, type=transaction_type
    ).save()

    yield currency_feature, date_feature, transaction_type, amount_usd, amount_eur

    # Cleanup
    amount_eur.delete()
    amount_usd.delete()
    transaction_type.delete()
    date_feature.delete()
    currency_feature.delete()


@pytest.fixture(scope="module")
def transactions_schema(transactions_features):
    # Create schema
    _, date_feature, _, amount_usd, amount_eur = transactions_features
    currency_feature = ln.Feature.get(name="currency_name")

    schema = ln.Schema(
        name="transaction_dataframe",
        otype="DataFrame",
        features=[
            date_feature,
            amount_usd,
            amount_eur,
            currency_feature,
        ],
    ).save()

    yield schema

    # Cleanup
    schema.delete()


@pytest.fixture
def transactions_dataframe():
    # Create sample data
    data = {
        "date": [
            datetime.date(2024, 1, 1),
            datetime.date(2024, 1, 2),
            datetime.date(2024, 1, 3),
            datetime.date(2024, 1, 4),
            datetime.date(2024, 1, 5),
        ],
        "transaction_amount_usd_cent": [1000, 2000, 3000, 4000, 5000],
        "transaction_amount_eur_cent": [850, 1700, 2550, 3400, 4250],
        "currency_name": ["USD", "EUR", "USD", "EUR", "USD"],
    }
    return pd.DataFrame(data)


def test_schema_creation(transactions_schema):
    """Test if schema was created properly"""
    schema = ln.Schema.get(name="transaction_dataframe")
    assert schema is not None
    assert schema.otype == "DataFrame"
    # check the order of the features
    assert schema.members.list("name") == [
        "date",
        "transaction_amount_usd_cent",
        "transaction_amount_eur_cent",
        "currency_name",
    ]


def test_currency_labels(currency_labels):
    """Test if currency labels were created properly"""
    currency_type, usd, eur = currency_labels

    # Test currency type
    assert ln.ULabel.get(name="Currency", is_type=True) is not None

    # Test currency labels
    assert ln.ULabel.get(name="USD") is not None
    assert ln.ULabel.get(name="EUR") is not None

    # Test relationships
    assert usd.type == currency_type
    assert eur.type == currency_type


def test_data_curation(transactions_schema, transactions_dataframe):
    """Test if data curation works properly"""
    curator = ln.curators.DataFrameCurator(transactions_dataframe, transactions_schema)
    assert curator.validate() is True
    artifact = curator.save_artifact(key="test_transaction_dataset.parquet")
    artifact.delete(permanent=True)


def test_missing_required_feature(transactions_schema):
    """Test if validation fails for invalid data"""
    data_missing_required_feature = {
        "date": [datetime.date(2024, 1, 1)],
        "transaction_amount_usd_cent": [1000],
        "currency_name": ["USD"],
    }
    invalid_df = pd.DataFrame(data_missing_required_feature)

    schema = ln.Schema.get(name="transaction_dataframe")
    curator = ln.curators.DataFrameCurator(invalid_df, schema)

    failure_code = curator.validate()
    assert (
        failure_code
        == "column 'transaction_amount_eur_cent' not in dataframe. Columns in dataframe: ['date', 'transaction_amount_usd_cent', 'currency_name']"
    )


def test_invalid_label(transactions_schema):
    """Test if validation fails for invalid currency"""
    # Create dataframe with invalid currency
    invalid_data = {
        "date": [datetime.date(2024, 1, 1)],
        "transaction_amount_usd_cent": [1000],
        "transaction_amount_eur_cent": [850],
        "currency_name": ["GBP"],  # Invalid currency not in our labels
    }
    invalid_df = pd.DataFrame(invalid_data)

    schema = ln.Schema.get(name="transaction_dataframe")
    curator = ln.curators.DataFrameCurator(invalid_df, schema)

    failure_code = curator.validate()
    assert "1 term is not validated: 'GBP'" in failure_code
