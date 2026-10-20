import lamindb as ln


def test_curator_remote():
    lamindata_artifacts = ln.Artifact.connect("laminlabs/lamindata")
    curator = ln.curators.DataFrameCurator(
        lamindata_artifacts.get("Ywz5JiVNHOWSJDiK"),
        schema=ln.examples.schemas.a_feature_type(),
    )
    curator.validate()
