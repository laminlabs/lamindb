from datetime import datetime

import pytest
from django.core.exceptions import ValidationError

# Assuming these models exist in your project
from lamindb.models import (
    Feature,
    FlexTable,
    Param,
    Run,
    RunData,
    Transform,
)

pytestmark = pytest.mark.django_db


@pytest.fixture
def feature():
    return Feature(name="test_feature", dtype="float").save()


@pytest.fixture
def param():
    return Param(name="test_param", dtype="str").save()


@pytest.fixture
def run():
    transform = Transform(key="test_transform").save()
    run = Run(name="Test Run", transform=transform).save()
    yield run
    run.delete()
    transform.delete()


@pytest.fixture
def tidy_table():
    return FlexTable(
        uid="test-table",
        name="Test Table",
        description="Test Description",
    ).save()


def test_run_data_int_value(run, feature):
    data = RunData(run=run, feature=feature, row=1, value_int=42).save()
    assert data.value_int == 42
    assert data.value_float is None


def test_run_data_float_value(run, feature):
    data = RunData(run=run, feature=feature, row=2, value_float=3.14).save()
    assert data.value_float == 3.14
    assert data.value_int is None


def test_run_data_string_value(run, feature):
    data = RunData(run=run, feature=feature, row=3, value_str="test string").save()
    assert data.value_str == "test string"


def test_run_data_datetime_value(run, feature):
    test_datetime = datetime(2025, 1, 1, 12, 0)
    data = RunData(run=run, feature=feature, row=4, value_datetime=test_datetime).save()
    assert data.value_datetime == test_datetime


def test_run_data_result_convention(run, feature):
    """Test that row=-1 can be used for result data"""
    result_data = RunData(run=run, feature=feature, row=-1, value_float=99.9).save()
    assert result_data.row == -1
    assert result_data.value_float == 99.9


def test_feature_param_mutual_exclusivity(run, feature, param):
    """Test that a row can't have both feature and param"""
    # Test with feature only
    feature_data = RunData(run=run, feature=feature, row=1, value_int=42).save()
    assert feature_data.feature == feature
    assert feature_data.param is None

    # Test with param only
    param_data = RunData(run=run, param=param, row=2, value_int=43).save()
    assert param_data.param == param
    assert param_data.feature is None

    # Test with both feature and param
    with pytest.raises(ValidationError):
        data = RunData(run=run, feature=feature, param=param, row=3, value_int=44)
        data.full_clean() if hasattr(data, "full_clean") else data._full__clean()
