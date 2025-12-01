import lamindb as ln
import pytest


def test_querydb_multiple_instances():
    """Accessing multiple instances simultaneously must work."""
    cxg = ln.QueryDB("laminlabs/cellxgene")
    lamindata = ln.QueryDB("laminlabs/lamindata")
    qs1 = cxg.artifacts.filter(suffix=".h5ad")
    qs2 = lamindata.artifacts.filter(suffix=".zarr")
    assert qs1._db != qs2._db


def test_querydb_bionty():
    """Querying a record from bionty must work."""
    cxg = ln.QueryDB("laminlabs/cellxgene")
    assert len(cxg.genes.filter(symbol__startswith="TP53")) > 0


def test_querydb_missing_module():
    """Attempting to access an attribute that comes from a missing module must error."""
    db = ln.QueryDB("laminlabs/lamin-site-assets")  # instance without bionty

    with pytest.raises(AttributeError) as e:
        db.genes.first()

    assert (
        "Registry 'genes' not found in installed modules for instance 'laminlabs/lamin-site-assets'."
        in str(e.value)
    )


def test_querydb_rejects_capitalized():
    """Accessing registries with capitalized names must fail."""
    cxg = ln.QueryDB("laminlabs/cellxgene")
    with pytest.raises(AttributeError) as e:
        cxg.Artifacts  # noqa: B018
    assert "Use lowercase plural form" in str(e.value)
