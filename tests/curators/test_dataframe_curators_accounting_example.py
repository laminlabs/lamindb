"""Test suite for accounting on bank transactions."""

import datetime

import lamindb as ln
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def transactions_schema():
    # Labels
    currency_type = ln.ULabel(name="Currency", is_type=True).save()
    usd = ln.ULabel(name="USD", type=currency_type).save()
    eur = ln.ULabel(name="EUR", type=currency_type).save()

    assert usd.type == currency_type
    assert eur.type == currency_type

    # Features
    currency = ln.Feature(name="currency_name", dtype="cat[ULabel[Currency]]").save()
    date = ln.Feature(name="date", dtype="date").save()

    transaction_type = ln.Feature(name="Transaction", is_type=True).save()
    amount_usd = ln.Feature(
        name="transaction_amount_usd_cent", dtype=int, type=transaction_type
    ).save()
    amount_eur = ln.Feature(
        name="transaction_amount_eur_cent", dtype=int, type=transaction_type
    ).save()

    # Schema
    schema = ln.Schema(
        name="transaction_dataframe",
        otype="DataFrame",
        features=[
            date,
            amount_usd,
            amount_eur,
            currency,
        ],
        coerce_dtype=True,
    ).save()

    yield schema

    schema.delete()
    amount_eur.delete()
    amount_usd.delete()
    transaction_type.delete()
    date.delete()
    currency.delete()
    eur.delete()
    usd.delete()
    currency_type.delete()


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


def test_data_curation(transactions_schema, transactions_dataframe):
    """Test if data curation works properly"""
    curator = ln.curators.DataFrameCurator(transactions_dataframe, transactions_schema)
    assert curator.validate() is None
    artifact = curator.save_artifact(key="test_transaction_dataset.csv")
    assert artifact.suffix == ".csv"
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

    with pytest.raises(ln.errors.ValidationError) as err:
        curator.validate()
        message = "column 'transaction_amount_eur_cent' not in dataframe. Columns in dataframe: ['date', 'transaction_amount_usd_cent', 'currency_name']"
        assert str(err) == message
        assert err.exconly() == f"lamindb.errors.ValidationError: {message}"


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

    with pytest.raises(ln.errors.ValidationError):
        curator.validate()
    # exconly = """lamindb.errors.ValidationError: 1 term is not validated: 'GBP'
    # → fix typos, remove non-existent values, or save terms via .add_new_from("currency_name")"""
    # assert err.exconly() == exconly
