import lamindb as ln
import pytest
from lamindb.errors import ValidationError


def test_feature_regex_transfer():
    handle = ln.setup.settings.user.handle
    ln.connect("testdb1")
    feature = ln.Feature(name="transfer_regex_sample", dtype=str, regex="A+").save()
    feature._aux["other_metadata"] = True
    feature.save()
    ln.connect("testdb2")
    ln.DB(f"{handle}/testdb1").Feature.get(uid=feature.uid).save()
    transferred = ln.Feature.get(uid=feature.uid)
    assert transferred.regex == "A+"
    assert transferred._aux["other_metadata"] is True
    from lamindb.curators.core import DictCurator

    schema = ln.Schema([transferred]).save()
    DictCurator({transferred.name: "AA"}, schema).validate()
    with pytest.raises(ValidationError, match="regex"):
        DictCurator({transferred.name: "AB"}, schema).validate()
