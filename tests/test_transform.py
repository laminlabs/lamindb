import lamindb as ln


def test_create():
    transform = ln.Transform(stem_id="NJvdsWWbJlZS", version="0")
    assert transform.id == "NJvdsWWbJlZSz8"
