import lamindb as ln
import pytest


def test_querydb_multiple_instances():
    """Accessing multiple instances simultaneously must work."""
    cxg_db = ln.QueryDB("laminlabs/cellxgene")
    lamindata_db = ln.QueryDB("laminlabs/lamindata")
    qs1 = cxg_db.Artifact.filter(suffix=".h5ad")
    qs2 = lamindata_db.Artifact.filter(suffix=".zarr")
    assert qs1._db != qs2._db


def test_querydb_bionty():
    """Querying a record from bionty must work."""
    cxg_db = ln.QueryDB("laminlabs/cellxgene")
    assert len(cxg_db.bionty.Gene.filter(symbol__startswith="TP53")) > 0


def test_querydb_missing_module():
    """Attempting to access an attribute that comes from a missing module must error."""
    site_assets_db = ln.QueryDB(
        "laminlabs/lamin-site-assets"
    )  # instance without bionty

    with pytest.raises(AttributeError) as e:
        site_assets_db.bionty.Gene.first()

    assert (
        "Schema 'bionty' not available in instance 'laminlabs/lamin-site-assets'."
        in str(e.value)
    )


def test_querydb_instantiate_class():
    """Attempting to instantiate a class must error."""
    cxg_db = ln.QueryDB("laminlabs/cellxgene")
    with pytest.raises(TypeError) as e:
        cxg_db.Artifact()
    assert (
        "Cannot instantiate Artifact from QueryDB. Use Artifact.filter(), Artifact.get(), etc. to query records."
        in str(e.value)
    )


def test_querydb_rejects_lowercase():
    """Accessing registries with lowercase names must fail."""
    cxg_db = ln.QueryDB("laminlabs/cellxgene")
    with pytest.raises(AttributeError) as e:
        cxg_db.artifacts  # noqa: B018
    assert "Registry names must be capitalized and singular." in str(e.value)


def test_querydb_cache():
    """Subsequent accesses must return cached wrapper."""
    cxg_db = ln.QueryDB("laminlabs/cellxgene")
    artifact1 = cxg_db.Artifact
    artifact2 = cxg_db.Artifact
    assert artifact1 is artifact2


def test_querydb_dir():
    """__dir__ must return discovered registries."""
    cxg = ln.QueryDB("laminlabs/cellxgene")
    dir_result = dir(cxg)
    assert "Artifact" in dir_result
    assert "Collection" in dir_result
    assert "Gene" not in dir_result
    assert "bt" in dir_result
