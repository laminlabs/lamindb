import lamindb as ln


def test_view():
    ln.view(modules="core")
    ln.view()
