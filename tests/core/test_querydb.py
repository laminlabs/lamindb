import lamindb as ln
import pytest


def test_DB_multiple_instances():
    """Accessing multiple instances simultaneously must work."""
    cxg_db = ln.DB("laminlabs/cellxgene")
    lamindata_db = ln.DB("laminlabs/lamindata")
    qs1 = cxg_db.Artifact.filter(suffix=".h5ad")
    qs2 = lamindata_db.Artifact.filter(suffix=".zarr")
    assert qs1._db != qs2._db


def test_DB_bionty():
    """Querying a record from bionty must work."""
    cxg_db = ln.DB("laminlabs/cellxgene")
    assert len(cxg_db.bionty.Gene.filter(symbol__startswith="TP53")) > 0


def test_DB_missing_module():
    """Attempting to access an attribute that comes from a missing module must error."""
    site_assets_db = ln.DB("laminlabs/lamin-site-assets")  # instance without bionty

    with pytest.raises(AttributeError) as e:
        site_assets_db.bionty.Gene.first()

    assert (
        "Schema 'bionty' not available in instance 'laminlabs/lamin-site-assets'."
        in str(e.value)
    )


def test_DB_instantiate_class():
    """Attempting to instantiate a class must error."""
    cxg_db = ln.DB("laminlabs/cellxgene")
    with pytest.raises(TypeError) as e:
        cxg_db.Artifact()
    assert (
        "Cannot instantiate Artifact from DB. Use Artifact.filter(), Artifact.get(), etc. to query records."
        in str(e.value)
    )


@pytest.mark.parametrize(
    "attr,expected_msg",
    [
        ("artifacts", "Registry 'artifacts' not found"),
        ("foo", "Registry 'foo' not found"),
        ("celltype", "Registry 'celltype' not found"),
    ],
)
def test_DB_rejects_invalid_attributes(attr, expected_msg):
    """Accessing invalid attributes must fail."""
    cxg_db = ln.DB("laminlabs/cellxgene")
    with pytest.raises(AttributeError) as e:
        getattr(cxg_db, attr)
    assert expected_msg in str(e.value)


def test_DB_cache():
    """Subsequent accesses must return cached wrapper."""
    cxg_db = ln.DB("laminlabs/cellxgene")
    artifact1 = cxg_db.Artifact
    artifact2 = cxg_db.Artifact
    assert artifact1 is artifact2


def test_queryset_caching():
    """Calling `.filter()` multiple times should return different results."""
    cxg_db = ln.DB("laminlabs/cellxgene")
    res_1 = cxg_db.Artifact.filter().first()
    res_2 = cxg_db.Artifact.filter().last()

    assert res_1 != res_2


def test_DB_dir():
    """__dir__ must return discovered registries."""
    cxg = ln.DB("laminlabs/cellxgene")
    dir_result = dir(cxg)
    assert "Artifact" in dir_result
    assert "Collection" in dir_result
    assert "Gene" not in dir_result
    assert "bionty" in dir_result
