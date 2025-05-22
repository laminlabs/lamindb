from unittest.mock import patch

import bionty as bt
import lamindb as ln
import pytest
from lamindb.models._django import get_artifact_with_related


def test_transfer_from_remote_to_local():
    """Test transfer from remote to local instance."""

    bt.Gene.filter().delete()
    bt.Organism.filter().delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()

    # test transfer form an instance with extra registry modules (laminlabs/lamin-dev)
    # transfer 1st artifact
    artifact = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(uid="livFRRpMaOgb3y8U2mK2")
        .one()
    )

    # test describe postgres
    result = get_artifact_with_related(
        artifact,
        include_m2m=True,
        include_fk=True,
        include_feature_link=True,
        include_schema=True,
    )
    assert result["related_data"]["m2m"]["tissues"] == {2: "cortex of kidney"}
    assert sorted(
        result["related_data"]["link"]["links_ulabel"], key=lambda d: d["id"]
    ) == [
        {"id": 7, "ulabel": 15, "feature": 1},
        {"id": 8, "ulabel": 10, "feature": 10},
    ]
    assert result["related_data"]["schemas"][615][0] == "obs"
    assert result["related_data"]["schemas"][615][1] == {
        "Feature": [
            "donor_id",
            "development_stage",
            "disease",
            "cell_type",
            "sex",
            "assay",
            "tissue",
            "self_reported_ethnicity",
            "tissue_type",
            "suspension_type",
            "organism",
        ]
    }
    assert result["related_data"]["fk"]["storage"] == {
        "id": 4,
        "name": "s3://cellxgene-data-public",
    }

    id_remote = artifact.id
    run_remote = artifact.run
    transform_remote = artifact.transform
    created_by_remote = artifact.created_by
    storage_remote = artifact.storage
    organism_remote = artifact.organisms.get(name="human")

    artifact.save()

    # check all ids are adjusted
    assert id_remote != artifact.id
    assert run_remote != artifact.run
    assert transform_remote != artifact.transform
    assert created_by_remote.handle != artifact.created_by.handle
    assert storage_remote.uid == artifact.storage.uid
    assert storage_remote.created_at == artifact.storage.created_at
    organism = artifact.organisms.get(name="human")
    assert organism.created_at != organism_remote.created_at

    # now check that this is idempotent and we can run it again
    artifact_repeat = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(uid="livFRRpMaOgb3y8U2mK2")
        .one()
    )
    artifact_repeat.save()

    # now prepare a new test case
    # mimic we have an existing feature with a different uid but same name
    feature = ln.Feature.get(name="organism")
    feature.uid = "existing"
    feature.save()

    # transfer 2nd artifact
    artifact2 = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(uid="qz35YaRk09XtYAyLvjyZ")
        .one()
    )
    artifact2.save()

    # check the feature name
    bt.settings.organism = "mouse"
    assert artifact2.organisms.get(name="mouse") == bt.settings.organism
    assert artifact.features["obs"].get(name="organism").uid == "existing"

    # test transfer form an instance with fewer modules (laminlabs/lamin-site-assets)
    artifact3 = ln.Artifact.using("laminlabs/lamin-site-assets").get(
        "lgRNHNtMxjU0y8nIagt7"
    )
    # load saves the artifact to default instance
    artifact3.load()

    ln.Artifact.filter().delete(permanent=True, storage=False)
    ln.Schema.filter().delete()
    bt.Gene.filter().delete()
    bt.Organism.filter().delete()
    ln.ULabel.filter().delete()
    bt.Disease.filter().delete()
    bt.CellLine.filter().delete()
    bt.CellType.filter().delete()
    bt.Phenotype.filter().delete()
    bt.Ethnicity.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.DevelopmentalStage.filter().delete()
    bt.Tissue.filter().delete()
    ln.Feature.filter().delete()


def test_transfer_into_space():
    # grab any ulabel from the default space
    ulabel = ln.ULabel.using("laminlabs/lamin-dev").filter(space__id=1).first()

    space = ln.Space(name="space for transfer", uid="00000123").save()
    with patch.object(ln.context, "_space", new=space):
        ulabel.save()
    assert ulabel.space_id == space.id

    ulabel.delete()
    space.delete()


def test_using_record_organism():
    """Test passing record and organism to the using_key instance."""
    import bionty as bt

    release_110 = bt.Source.filter(
        organism="mouse", entity="bionty.Gene", version="release-110"
    ).one()
    release_110_cxg = (
        bt.Source.using("laminlabs/lamin-dev")
        .filter(organism="mouse", entity="bionty.Gene", version="release-110")
        .one()
    )
    release_112_cxg = (
        bt.Source.using("laminlabs/lamin-dev")
        .filter(organism="mouse", entity="bionty.Gene", version="release-112")
        .one()
    )

    # passing the wrong source
    inspector = bt.Gene.using("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_112_cxg,
        strict_source=True,
    )
    assert len(inspector.validated) == 0

    # passing the correct source
    inspector = bt.Gene.using("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_110_cxg,
        strict_source=True,
    )
    assert len(inspector.validated) == 2

    # passing the correct source but from the wrong instance
    with pytest.raises(ValueError) as error:
        inspector = bt.Gene.using("laminlabs/lamin-dev").inspect(
            ["ENSMUSG00000102862", "ENSMUSG00000084826"],
            field=bt.Gene.ensembl_gene_id,
            source=release_110,
        )
    assert (
        "record must be a bionty.Source record from instance 'laminlabs/lamin-dev'"
        in str(error.value)
    )
