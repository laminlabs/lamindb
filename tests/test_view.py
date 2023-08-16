import lamindb as ln


def test_vew():
    ln.view(schema="core")
    ln.view()
