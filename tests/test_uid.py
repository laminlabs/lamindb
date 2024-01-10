import lamindb as ln
import lnschema_bionty as lb
import pytest
from lnschema_bionty._bionty import encode_uid


def test_lb_encode_uid():
    assert (
        encode_uid(lb.Gene, {"ensembl_gene_id": "ENSG00000081059", "symbol": "TCF7"})[
            "uid"
        ]
        == "7IkHKPl0ScQR"
    )
    with pytest.raises(AssertionError):
        encode_uid(lb.Organism, {"ensembl_gene_id": "ENSG00000081059"})
    assert encode_uid(lb.CellType, {"ontology_id": "CL:0000084"})["uid"] == "22LvKd01"
    assert (
        encode_uid(lb.Organism, {"ontology_id": "NCBITaxon:9606", "name": "human"})[
            "uid"
        ]
        == "1dpCL6Td"
    )
    assert encode_uid(lb.Organism, {"name": "human"})["uid"] == "4gQdjtxb"
