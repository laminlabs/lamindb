import lamindb as ln
import pytest


@pytest.mark.parametrize("host_type", ["artifact", "run"])
@pytest.mark.parametrize(
    "value", ["first part, second part", 'a "quoted", value\nnext line', ""]
)
def test_text_feature_roundtrip(host_type, value):
    feature = ln.Feature(name="text_roundtrip_description", dtype=str).save()
    transform = None
    if host_type == "artifact":
        host = ln.Artifact(".gitignore", key="text_roundtrip").save()
    else:
        transform = ln.Transform(key="text_roundtrip").save()
        host = ln.Run(transform).save()
    try:
        host.features.add_values({feature.name: value})
        assert host.features.get_values() == {feature.name: value}
    finally:
        host.delete(permanent=True)
        if transform is not None:
            transform.delete(permanent=True)
        feature.delete(permanent=True)
