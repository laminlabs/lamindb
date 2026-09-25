"""Transfer from laminlabs/lamin-dev. That instance is private and needs a test user."""

from unittest.mock import patch

import bionty as bt
import lamindb as ln
import pytest
from lamindb.models._django import get_artifact_or_run_with_related


@pytest.fixture(autouse=True)
def _require_testuser():
    handle = ln.setup.settings.user.handle
    if handle not in {"testuser1", "testuser2"}:
        pytest.skip("laminlabs/lamin-dev requires testuser authentication")


def test_transfer_from_lamin_dev(connected_bionty):
    """Test transfer from the private lamin-dev instance."""
    # lamin-dev has an extra schema module, pertdb, and this artifact is labeled with it.
    artifact1 = ln.Artifact.connect("laminlabs/lamin-dev").get("livFRRpM")

    result = get_artifact_or_run_with_related(
        artifact1,
        include_m2m=True,
        include_fk=True,
        include_feature_link=True,
        include_schema=True,
    )
    assert result["related_data"]["m2m"]["tissues"] == {
        2: {
            "id": 2,
            "uid": "6VHBo6XsJZqmaQ",
            "abbr": None,
            "name": "cortex of kidney",
            "tissue": 2,
            "feature": None,
            "ontology_id": "UBERON:0001225",
            "tissue_display": "cortex of kidney",
        }
    }
    assert sorted(
        result["related_data"]["link"]["links_ulabel"], key=lambda d: d["id"]
    ) == [
        {
            "id": 7,
            "uid": "ydyPUMjh",
            "name": "donor_24",
            "ulabel": 15,
            "feature": 1,
            "reference": None,
            "reference_type": None,
            "ulabel_display": "donor_24",
        },
        {
            "id": 8,
            "uid": "JJ3d8a2v",
            "name": "na",
            "ulabel": 10,
            "feature": 10,
            "reference": None,
            "reference_type": None,
            "ulabel_display": "na",
        },
    ]
    assert result["related_data"]["m2m_schemas"][615][0] == "obs"
    assert result["related_data"]["m2m_schemas"][615][1] == {
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

    # slot "tada" is a wetlab.Compound schema ("test schema for triggering error in CI")
    with pytest.raises(ValueError, match="schema slot 'tada'"):
        artifact1.save(transfer="annotations")
    if artifact1.pk is not None and ln.Artifact.filter(uid=artifact1.uid).exists():
        artifact1.delete(storage=False, permanent=True)

    # An artifact whose schemas the test instance can load.
    artifact2 = ln.Artifact.connect("laminlabs/lamin-dev").get("qz35YaRk")
    id_remote = artifact2.id
    run_remote = artifact2.run
    transform_remote = artifact2.transform
    created_by_remote = artifact2.created_by
    storage_remote = artifact2.storage
    organism_remote = artifact2.organisms.get(name="mouse")

    artifact2.save(transfer="annotations")

    assert id_remote != artifact2.id
    if run_remote is not None:
        assert artifact2.run is not None
        assert run_remote.uid != artifact2.run.uid
    if transform_remote is not None:
        assert artifact2.transform is not None
        assert transform_remote.uid != artifact2.transform.uid
    assert created_by_remote.uid == artifact2.created_by.uid
    assert created_by_remote.handle == artifact2.created_by.handle
    assert storage_remote.uid == artifact2.storage.uid
    assert storage_remote.created_at == artifact2.storage.created_at
    organism = artifact2.organisms.get(name="mouse")
    assert organism.uid == organism_remote.uid
    assert organism._state.db in {None, "default"}

    artifact_repeat = ln.Artifact.connect("laminlabs/lamin-dev").get(artifact2.uid)
    artifact_repeat.save(transfer="annotations")
    assert artifact_repeat.id == artifact2.id

    # An existing feature with the same name keeps its uid across a re-transfer.
    feature = artifact2.features.slots["obs"].members.get(name="tissue")
    existing_uid = f"exst{feature.uid[-8:]}"
    feature.uid = existing_uid
    feature.save()
    ln.Artifact.connect("laminlabs/lamin-dev").get(artifact2.uid).save(
        transfer="annotations"
    )
    assert (
        artifact2.features.slots["obs"].members.get(name="tissue").uid == existing_uid
    )

    artifact2.delete(storage=False, permanent=True)


def test_transfer_keeps_source_space(connected_bionty):
    # A source object in the shared `all` space stays there. The current space
    # must not replace it.
    ulabel = (
        ln.ULabel.connect("laminlabs/lamin-dev").filter(space__uid="a" * 12).first()
    )
    source_space_uid = ulabel.space.uid

    space = ln.Space(name="space for transfer", uid="00000123").save()
    with patch.object(ln.context, "_space", new=space):
        ulabel.save()
    assert ulabel.space.uid == source_space_uid
    assert ulabel.space_id != space.id

    ulabel.delete(permanent=True)
    ln.Run.filter(space=space).delete(permanent=True)
    ln.Transform.filter(space=space).delete(permanent=True)
    space.delete()


def test_using_record_organism(connected_bionty):
    """Test passing record and organism to the using_key instance."""
    release_110_cxg = bt.Source.connect("laminlabs/lamin-dev").get(
        organism="mouse", entity="bionty.Gene", version="release-110"
    )
    release_112_cxg = bt.Source.connect("laminlabs/lamin-dev").get(
        organism="mouse", entity="bionty.Gene", version="release-112"
    )
    release_110 = release_110_cxg.save()  # transfer source record
    release_110_cxg = bt.Source.connect("laminlabs/lamin-dev").get(
        organism="mouse", entity="bionty.Gene", version="release-110"
    )

    inspector = bt.Gene.connect("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_112_cxg,
        strict_source=True,
    )
    assert len(inspector.validated) == 0

    inspector = bt.Gene.connect("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_110_cxg,
        strict_source=True,
    )
    assert len(inspector.validated) == 2

    with pytest.raises(ValueError) as error:
        bt.Gene.connect("laminlabs/lamin-dev").inspect(
            ["ENSMUSG00000102862", "ENSMUSG00000084826"],
            field=bt.Gene.ensembl_gene_id,
            source=release_110,
        )
    assert (
        "record must be a bionty.Source record from instance 'laminlabs/lamin-dev'"
        in str(error.value)
    )
