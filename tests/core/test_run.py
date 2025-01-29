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
    run = ln.Run(transform)
    assert run.reference is None
    assert run.reference_type is None
    run2 = ln.Run(transform, reference="test1", reference_type="test2")
    assert run2.reference == "test1"
    assert run2.reference_type == "test2"
    assert run.uid != run2.uid
    transform.delete()


def test_edge_cases():
    with pytest.raises(ValueError) as error:
        ln.Run(1, 2)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: transform"
    with pytest.raises(TypeError) as error:
        ln.Run()
    assert error.exconly() == "TypeError: Pass transform parameter"
