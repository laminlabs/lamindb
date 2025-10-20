import lamindb as ln


def test_curator_remote():
    lamindata_artifacts = ln.Artifact.using("laminlabs/lamindata")
    curator = ln.curators.DataFrameCurator(
        lamindata_artifacts.get("F6PrJLEnWeHqUHrK0000"),
        schema=ln.examples.schemas.valid_features(),
    )
    curator.validate()
