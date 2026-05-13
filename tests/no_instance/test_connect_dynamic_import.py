def test_connect_dynamic_import(ccaplog):
    import lamindb as ln

    # this only currently works if not instance was configured in the environment
    # in all other cases, we still trigger a reset_django() and hence django variables
    # become stale in case of a dynamic import
    assert ln.setup.settings.instance.slug == "none/none"

    ln.connect("laminlabs/lamin-site-assets")
    assert "connected in read-only mode" in ccaplog.text
    assert ln.Artifact.filter(key__startswith="blog").count() > 0
    ln.setup.disconnect()
