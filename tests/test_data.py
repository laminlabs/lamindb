import lamindb as ln
from lamindb.dev._data import add_transform_to_kwargs


def test_add_transform_to_kwargs():
    kwargs = {}
    transform = ln.Transform(name="hello")
    run = ln.Run(transform=transform)
    add_transform_to_kwargs(kwargs, run)
    assert kwargs["transform"] == transform
