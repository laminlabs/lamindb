# also see the transfer.md guide (golden artifact + record) and tests/transfer

from unittest.mock import patch

import bionty as bt
import lamindb as ln
import pytest
from lamindb.models._django import get_artifact_or_run_with_related


def test_schema_transfer_defaults_to_annotations():
    db = ln.DB("laminlabs/lamindata")
    schema_uid = "pnQvQVcQ417bfmVq"
    remote_schema = db.Schema.get(schema_uid)
    remote_member_names = remote_schema.members.to_list("name")
    assert len(remote_member_names) > 0
    remote_perturbation_type_name = "Perturbation"

    existing_local = ln.Schema.filter(uid=schema_uid).one_or_none()
    if existing_local is not None:
        existing_local.delete(permanent=True)

    existing_local_perturbation_type = ln.Record.filter(
        name=remote_perturbation_type_name, is_type=True
    ).one_or_none()
    if existing_local_perturbation_type is not None:
        ln.Record.filter(type=existing_local_perturbation_type).delete(permanent=True)
        existing_local_perturbation_type.delete(permanent=True)

    assert (
        ln.Record.filter(name=remote_perturbation_type_name, is_type=True).count() == 0
    )
    assert ln.Record.filter(type__name=remote_perturbation_type_name).count() == 0

    record_only = db.Schema.get(schema_uid).save(transfer="record")
    assert record_only.members.count() == 0

    transferred = db.Schema.get(schema_uid).save()
    transferred_member_names = transferred.members.to_list("name")
    assert transferred_member_names == remote_member_names

    perturbation_feature = transferred.members.get(name="perturbation")
    assert perturbation_feature.dtype_as_str.startswith("cat[Record[")
    assert perturbation_feature.dtype_as_object is not None

    # The dtype only needs the Perturbation type. Its records stay behind.
    assert (
        ln.Record.filter(name=remote_perturbation_type_name, is_type=True).count() == 1
    )
    assert ln.Record.filter(type__name=remote_perturbation_type_name).count() == 0

    before_count = transferred.links_feature.count()
    transferred_repeat = db.Schema.get(schema_uid).save()
    assert transferred_repeat.id == transferred.id
    assert transferred_repeat.links_feature.count() == before_count


def test_describe_artifact_from_remote_instance(capsys):
    # test describing from a remote instance with less modules
    artifact = ln.Artifact.connect("laminlabs/lamin-site-assets").first()
    artifact.describe()
    captured = capsys.readouterr()
    assert len(captured.out) > 50
    assert "artifact" in captured.out.lower()


def test_transfer_from_remote_to_local(ccaplog):
    """Test transfer from remote to local instance."""

    # test transfer from an instance with an extra schema module: pertdb
    # we also made sure that the artifact here has a pertdb label attached

    # transfer 1st artifact
    artifact1 = ln.Artifact.connect("laminlabs/lamin-dev").get("livFRRpM")

    # test describe postgres
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

    # transfer an artifact whose schemas the test instance can load
    artifact2 = ln.Artifact.connect("laminlabs/lamin-dev").get("qz35YaRk")
    id_remote = artifact2.id
    run_remote = artifact2.run
    transform_remote = artifact2.transform
    created_by_remote = artifact2.created_by
    storage_remote = artifact2.storage
    organism_remote = artifact2.organisms.get(name="mouse")

    artifact2.save(transfer="annotations")

    # check all ids are adjusted
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

    # now check that this is idempotent and we can run it again
    artifact_repeat = ln.Artifact.connect("laminlabs/lamin-dev").get(artifact2.uid)
    artifact_repeat.save(transfer="annotations")
    assert artifact_repeat.id == artifact2.id

    # an existing feature with the same name keeps its uid across a re-transfer
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

    # test transfer from an instance with fewer modules (laminlabs/lamin-site-assets)
    artifact3 = ln.Artifact.connect("laminlabs/lamin-site-assets").get(
        "lgRNHNtMxjU0y8nIagt7"
    )
    # test that implicit saving through `load()` works (also occurs for `cache()` or `open()` for run input tracking)
    artifact3.load()

    # delete with storage=False, because these are all stored in the source instances
    artifact2.delete(storage=False, permanent=True)
    artifact3.delete(
        storage=False
    )  # there is an issue here with permanent deletion because of schema module mismatch


def test_transfer_keeps_source_space():
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


def test_transfer_missing_space_errors():
    from lamindb.errors import NoWriteAccess
    from lamindb.models.sqlrecord import update_fk_to_default_db

    missing = ln.Space(name="restricted-perturbations", uid="noattach1")
    missing.id = 99
    record = ln.Record(name="space gate")
    record.uid = "recSpace"
    record.space = missing

    with pytest.raises(NoWriteAccess, match="restricted-perturbations") as error:
        update_fk_to_default_db(
            record,
            "space",
            None,
            {"mapped": [], "transferred": [], "run": True},
        )
    message = str(error.value)
    target = ln.setup.settings.instance.slug
    assert (
        f"attach space 'restricted-perturbations' to the target database '{target}'"
        in message
    )
    assert "Record(uid='recSpace')" in str(error.value)


def test_using_record_organism():
    """Test passing record and organism to the using_key instance."""

    release_110_cxg = bt.Source.connect("laminlabs/lamin-dev").get(
        organism="mouse", entity="bionty.Gene", version="release-110"
    )
    release_112_cxg = bt.Source.connect("laminlabs/lamin-dev").get(
        organism="mouse", entity="bionty.Gene", version="release-112"
    )
    release_110 = release_110_cxg.save()  # transfer source record
    release_110_cxg = (  # re-fetch
        bt.Source.connect("laminlabs/lamin-dev").get(
            organism="mouse", entity="bionty.Gene", version="release-110"
        )
    )

    # passing the wrong source
    inspector = bt.Gene.connect("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_112_cxg,
        strict_source=True,
    )
    assert len(inspector.validated) == 0

    # passing the correct source
    inspector = bt.Gene.connect("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_110_cxg,
        strict_source=True,
    )
    assert len(inspector.validated) == 2

    # passing the correct source but from the wrong instance
    with pytest.raises(ValueError) as error:
        inspector = bt.Gene.connect("laminlabs/lamin-dev").inspect(
            ["ENSMUSG00000102862", "ENSMUSG00000084826"],
            field=bt.Gene.ensembl_gene_id,
            source=release_110,
        )
    assert (
        "record must be a bionty.Source record from instance 'laminlabs/lamin-dev'"
        in str(error.value)
    )


def test_annotation_transfer_requires_schema_module(monkeypatch):
    import lamindb_setup as ln_setup
    import pandas as pd

    feature = ln.Feature(
        name="organism_module_gate", dtype="cat[bionty.Organism]"
    ).save()
    record = ln.Record(name="module gate record").save()
    artifact = ln.Artifact.from_dataframe(
        pd.DataFrame({"a": [1]}), key="module-gate.parquet", description="module gate"
    ).save()
    schema = ln.Schema(name="organism module gate", itype=bt.Organism).save()
    artifact.schemas.add(schema, through_defaults={"slot": "var"})

    from lamindb.models.record import RecordJson

    RecordJson(record=record, feature=feature, value="human").save()

    instance = ln_setup.settings.instance
    modules = [module for module in instance.modules if module != "bionty"]
    monkeypatch.setattr(instance, "_schema_str", ",".join(modules))

    from lamindb.models.sqlrecord import transfer_record_feature_values

    try:
        with pytest.raises(ValueError, match="schema module"):
            transfer_record_feature_values(
                record,
                "default",
                record.pk,
                None,
                {"mapped": [], "transferred": [], "run": True},
            )
        with pytest.raises(ValueError, match="sqlrecord"):
            artifact.features._add_from(
                artifact, transfer_logs={"mapped": [], "transferred": [], "run": True}
            )
    finally:
        artifact.schemas.clear()
        artifact.delete(permanent=True)
        schema.delete(permanent=True)
        record.delete(permanent=True)
        feature.delete(permanent=True)


def test_using_query_by_feature():
    assert ln.Artifact.connect("laminlabs/cellxgene").filter(n_of_donors__gte=100)


def _source_user(uid: str, handle: str, name: str):
    user = ln.User(uid=uid, handle=handle, name=name)
    user._state.db = "laminlabs/lamindata"
    return user


def _transfer_logs():
    # A non-None run skips creating the transfer run in these unit tests.
    return {"mapped": [], "transferred": [], "run": True}


def test_map_user_annotation_uses_same_uid():
    from types import SimpleNamespace

    from lamindb.models.sqlrecord import _map_user_annotation

    uid = "usrAnnot"
    handle = "annot-user"
    existing = ln.User.filter(uid=uid).one_or_none()
    if existing is not None:
        existing.delete(permanent=True)

    source = _source_user(uid, handle, "Annot User")
    logs = _transfer_logs()
    try:
        assert _map_user_annotation(source, feature=None, transfer_logs=logs) == handle
        assert ln.User.filter(uid=uid).one().handle == handle
        assert ln.User.get(uid=uid).id != ln.setup.settings.user.id
        # a second sync maps the existing registry row
        assert _map_user_annotation(source, feature=None, transfer_logs=logs) == handle
        assert ln.User.filter(uid=uid).count() == 1
        feature = SimpleNamespace(_dtype_str="cat[User.uid]")
        assert _map_user_annotation(source, feature=feature, transfer_logs=logs) == uid
    finally:
        saved = ln.User.filter(uid=uid).one_or_none()
        if saved is not None:
            saved.delete(permanent=True)


def test_map_user_annotation_asks_to_add_collaborator(monkeypatch):
    from django.db import ProgrammingError
    from lamindb.errors import NoWriteAccess
    from lamindb.models.sqlrecord import _map_user_annotation

    uid = "usrDeny1"
    source = _source_user(uid, "denied-user", "Denied")

    def deny(self, *args, **kwargs):
        raise ProgrammingError(
            'new row violates row-level security policy for table "lamindb_user"'
        )

    monkeypatch.setattr(ln.User, "save", deny)
    with pytest.raises(NoWriteAccess, match="collaborator"):
        _map_user_annotation(source, feature=None, transfer_logs=_transfer_logs())
    assert ln.User.filter(uid=uid).one_or_none() is None
