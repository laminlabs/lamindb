import pandas as pd
from lamindb.models.query_manager import Lookup


def test_lookup():
    df = pd.DataFrame(
        {
            "name": ["Sample 1", "Sample 1", "sample 1", "1 sample", "", None],
            "meta1": [
                "metadata~1",
                "metadata~1~1",
                "metadata~1~1~1",
                "1 metadata",
                None,
                None,
            ],
        }
    )

    inst = Lookup(
        df=df, field="name", tuple_name="TestTuple", prefix="prefix", keep=False
    )
    lookup = inst.lookup()

    assert len(lookup.sample_1) == 3
    lookup_sample_1_dicts = [i._asdict() for i in lookup.sample_1]
    assert {"name": "Sample 1", "meta1": "metadata~1~1"} in lookup_sample_1_dicts
    assert {"name": "Sample 1", "meta1": "metadata~1"} in lookup_sample_1_dicts
    assert {"name": "sample 1", "meta1": "metadata~1~1~1"} in lookup_sample_1_dicts

    assert lookup.prefix_1_sample._asdict() == {
        "name": "1 sample",
        "meta1": "1 metadata",
    }

    lookup_dict = lookup.dict()
    assert len(lookup_dict) == 3
    assert isinstance(lookup_dict["Sample 1"], list)
    assert len(lookup_dict["Sample 1"]) == 2
    assert isinstance(lookup_dict["sample 1"], tuple)
    assert lookup_dict["1 sample"]._asdict() == {
        "name": "1 sample",
        "meta1": "1 metadata",
    }

    assert lookup.__class__.__name__ == "Lookup"
    assert lookup.prefix_1_sample.__class__.__name__ == "TestTuple"


def test_lookup_empty_df():
    assert Lookup(df=pd.DataFrame(columns=["name"]), field="name").lookup().dict() == {}


def test_lookup_multiple_records():
    df = pd.DataFrame(
        {
            "name": ["experiment", "experiment"],
            "value": ["value1", "value2"],
        }
    )

    lookup = Lookup(df=df, field="name", keep="first").lookup()
    assert lookup.experiment.value == "value1"
    assert type(lookup.experiment).__name__ == "MyTuple"

    lookup = Lookup(df=df, field="name", keep="last").lookup()
    assert lookup.experiment.value == "value2"
    assert type(lookup.experiment).__name__ == "MyTuple"

    lookup = Lookup(df=df, field="name", keep=False).lookup()
    assert len(lookup.experiment) == 2
    assert type(lookup.experiment).__name__ == "_ListValueWrapper"

    lookup = Lookup(df=df, field="name", keep="first").lookup(return_field="value")
    assert lookup.experiment == "value1"
    assert type(lookup.experiment) is str


def test_lookup_keyword_field():
    df = pd.DataFrame(
        {
            "name": ["DEL", "class", "normal"],
            "value": [1, 2, 3],
        }
    )

    lookup = Lookup(df=df, field="name").lookup()
    assert lookup.del_.value == 1
    assert lookup.class_.value == 2
    assert lookup.normal.value == 3
