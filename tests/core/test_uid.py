import bionty as bt
from bionty._bionty import encode_uid


def test_bionty_encode_uid():
    assert (
        encode_uid(bt.Gene, {"ensembl_gene_id": "ENSG00000081059", "symbol": "TCF7"})[
            "uid"
        ]
        == "7IkHKPl0ScQR"
    )
    assert encode_uid(bt.CellType, {"ontology_id": "CL:0000084"})["uid"] == "22LvKd01"
    assert (
        encode_uid(bt.Organism, {"ontology_id": "NCBITaxon:9606", "name": "human"})[
            "uid"
        ]
        == "1dpCL6Td"
    )
    assert encode_uid(bt.Organism, {"name": "human"})["uid"] == "4gQdjtxb"
    assert (
        encode_uid(
            bt.Source,
            {
                "entity": "Source",
                "name": "ensembl",
                "version": "release-112",
                "organism": "vertebrates",
            },
        )["uid"]
        == "5MUNZB0a"
    )
    bt.settings.organism = "human"
    assert (
        encode_uid(bt.CellMarker, {"name": "test", "organism": bt.settings.organism})[
            "uid"
        ]
        == "2dZ52W9noUDK"
    )
    bt.Gene.filter().all().delete()
    bt.settings.organism.delete()
