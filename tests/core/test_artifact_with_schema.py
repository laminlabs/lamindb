import lamindb as ln
import pytest


@pytest.mark.parametrize(
    "schema",
    [
        ln.Schema(
            features=[
                ln.Feature(name="species", dtype=str).save(),
                ln.Feature(name="split", dtype=str).save(),
            ]
        ),
        None,
    ],
)
def test_create_external_schema(tsv_file, schema):
    if schema:
        schema.save()
    else:
        (ln.Feature(name="split", dtype=str).save(),)
        (ln.Feature(name="species", dtype=str).save(),)
    artifact = ln.Artifact(
        tsv_file,
        features={"species": "bird", "split": "train"},
        schema=schema,
        description="test",
    ).save()
    assert artifact.features.get_values() == {"species": "bird", "split": "train"}

    artifact.delete(permanent=True)
    if schema:
        schema.delete(permanent=True)
    ln.Feature.get(name="species").delete(permanent=True)
    ln.Feature.get(name="split").delete(permanent=True)


def test_from_dataframe_external_schema(df):
    species = ln.Feature(name="species", dtype="str").save()
    split = ln.Feature(name="split", dtype="str").save()
    external_schema = ln.Schema(itype=ln.Feature).save()

    feat1 = ln.Feature(name="feat1", dtype="int").save()
    feat2 = ln.Feature(name="feat2", dtype="int").save()
    schema = ln.Schema(
        features=[feat1, feat2],
        slots={"__external__": external_schema},
        otype="DataFrame",
    ).save()

    artifact = ln.Artifact.from_dataframe(
        df,
        features={"species": "bird", "split": "train"},
        schema=schema,
        description="test dataframe with external features",
    ).save()

    assert artifact.features.get_values() == {"species": "bird", "split": "train"}

    # Cleanup
    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    external_schema.delete(permanent=True)

    for feature in [species, split, feat1, feat2]:
        feature.delete(permanent=True)
