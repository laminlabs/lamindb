"""Django proxy models of lamindb registries behave like the registry they proxy."""

import lamindb as ln
import pytest
from django.db import connection
from django.test.utils import CaptureQueriesContext
from lamindb.base.utils import concrete_model, get_registry_name
from lamindb.errors import ValidationError
from lamindb.models.sqlrecord import validate_fields


class ProxyArtifact(ln.Artifact):
    class Meta:
        proxy = True
        app_label = "lamindb_tests"


class ProxyProject(ln.Project):
    class Meta:
        proxy = True
        app_label = "lamindb_tests"


class ProxyULabel(ln.ULabel):
    class Meta:
        proxy = True
        app_label = "lamindb_tests"


class ProxySchema(ln.Schema):
    class Meta:
        proxy = True
        app_label = "lamindb_tests"


@pytest.fixture
def artifact():
    artifact = ln.Artifact("README.md", key="proxy-models/README.md").save()
    yield artifact
    artifact.delete(permanent=True)


def test_registry_name_resolves_proxies():
    assert concrete_model(ProxyArtifact) is ln.Artifact
    assert get_registry_name(ProxyArtifact) == "Artifact"
    assert get_registry_name(ln.Artifact) == "Artifact"
    assert ProxyArtifact.__get_name_with_module__() == "Artifact"
    assert ProxyArtifact.__get_module_name__() == "core"


def test_proxy_queries_validate_like_the_concrete_registry(artifact):
    assert "transform" in ProxyArtifact.__get_available_fields__()
    proxy = ProxyArtifact.filter(key=artifact.key).one()
    assert type(proxy) is ProxyArtifact
    assert proxy == artifact


def test_proxy_artifact_is_tracked_as_an_artifact_input(artifact):
    # Not ln.track(): under pytest it would claim the transform later test files track.
    transform = ln.Transform(key="proxy-models-input-tracking").save()
    run = ln.Run(transform).save()
    ln.context._run = run
    try:
        proxy = ProxyArtifact.get(artifact.id)
        proxy.cache(is_run_input=True)
        assert artifact in run.input_artifacts.all()
        assert not run.input_collections.exists()
    finally:
        ln.context._run = None
        run.delete(permanent=True)
        transform.delete(permanent=True)


def test_proxy_artifact_annotates_like_an_artifact(artifact):
    feature = ln.Feature(name="proxy_models_note", dtype=str).save()
    try:
        proxy = ProxyArtifact.get(artifact.id)
        proxy.features.add_values({"proxy_models_note": "hello"})
        assert artifact.features.get_values()["proxy_models_note"] == "hello"
        proxy.features.remove_values("proxy_models_note")
        assert "proxy_models_note" not in artifact.features.get_values()
    finally:
        feature.delete(permanent=True)


def test_deleting_a_proxy_type_deletes_its_children():
    parent = ProxyProject(name="proxy-models-type", is_type=True).save()
    child = ln.Project(name="proxy-models-child", type=parent).save()
    try:
        parent.delete()
        child.refresh_from_db()
        assert child.branch_id == -1
    finally:
        child.delete(permanent=True)
        parent.delete(permanent=True)


def test_proxy_uid_length_is_validated():
    label = ProxyULabel.__new__(ProxyULabel)
    with pytest.raises(ValidationError, match="`uid` must be exactly"):
        validate_fields(label, {"uid": "short", "name": "proxy-models"})
    # legacy 20-character Schema uids stay valid through a proxy
    schema = ProxySchema.__new__(ProxySchema)
    validate_fields(schema, {"uid": "0" * 20})
    with pytest.raises(ValidationError, match="must be exactly 16 characters"):
        validate_fields(schema, {"uid": "0" * 12})


def test_bulk_save_batches_proxies_with_their_registry():
    labels = [ProxyULabel(name="proxy-models-a"), ln.ULabel(name="proxy-models-b")]
    try:
        with CaptureQueriesContext(connection) as queries:
            ln.save(labels)
        inserts = [
            q for q in queries if q["sql"].startswith('INSERT INTO "lamindb_ulabel"')
        ]
        assert len(inserts) == 1
        assert ln.ULabel.filter(name__startswith="proxy-models-").count() == 2
    finally:
        ln.ULabel.filter(name__startswith="proxy-models-").delete(permanent=True)


def test_deleting_a_proxy_artifact_removes_row_and_file(tmp_path):
    source = tmp_path / "proxy-delete.txt"
    source.write_text("proxy-models delete")
    uid = ln.Artifact(source, key="proxy-models/delete.txt").save().uid
    proxy = ProxyArtifact.get(uid)
    path = proxy.path
    assert path.exists()
    proxy.delete(permanent=True, storage=True)
    assert not ln.Artifact.filter(uid=uid).exists()
    assert not path.exists()
