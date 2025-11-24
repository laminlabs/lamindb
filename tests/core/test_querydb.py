import lamindb as ln


def test_querydb_multiple_instances():
    """Test accessing multiple instances simultaneously."""
    cellxgene = ln.QueryDB("laminlabs/cellxgene")
    lamindata = ln.QueryDB("laminlabs/lamindata")
    qs1 = cellxgene.artifacts.filter(suffix=".h5ad")
    qs2 = lamindata.artifacts.filter(suffix=".zarr")
    assert qs1._db != qs2._db


def test_querydb_bionty():
    """Test querying a record that is not from LaminDB like bionty."""
    cxg = ln.QueryDB("laminlabs/cellxgene")
    assert len(cxg.genes.filter(symbol__startswith="TP53")) > 0
