import lamindb as ln
from lamindb.dev._data import add_transform_to_kwargs


def test_add_transform_to_kwargs():
    kwargs = {}
    transform = ln.Transform(name="hello")
    transform.save()
    run = ln.Run(transform)
    add_transform_to_kwargs(kwargs, run)
    assert kwargs["transform"] == transform
