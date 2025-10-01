import lamindb as ln
import pytest


def test_run():
    transform = ln.Transform(key="My transform")
    with pytest.raises(ValueError) as error:
        ln.Run(transform)
    assert (
        error.exconly()
        == "ValueError: Please save transform record before creating a run"
    )
    transform.save()
    run = ln.Run(transform).save()
    assert run.status == "scheduled"
    assert run.reference is None
    assert run.reference_type is None
    run2 = ln.Run(transform, reference="test1", reference_type="test2").save()
    assert run2.reference == "test1"
    assert run2.reference_type == "test2"
    assert run.uid != run2.uid

    report_artifact = ln.Artifact("README.md", description="report of run2").save()
    run2.report = report_artifact
    environment = ln.Artifact("CONTRIBUTING.md", description="env of run2").save()
    run2.environment = environment

    run2.delete(permanent=True)

    # test deletion of run including attached artifacts
    assert ln.Artifact.objects.filter(uid=report_artifact.uid).exists() is False
    assert ln.Artifact.objects.filter(uid=environment.uid).exists() is False

    transform.delete(permanent=True)

    assert ln.Run.filter(uid=run.uid).count() == 0


def test_edge_cases():
    with pytest.raises(ValueError) as error:
        ln.Run(1, 2)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: transform"
    with pytest.raises(TypeError) as error:
        ln.Run()
    assert error.exconly() == "TypeError: Pass transform parameter"
