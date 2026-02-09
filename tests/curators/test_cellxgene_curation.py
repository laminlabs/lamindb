from typing import Generator

import bionty as bt
import lamindb as ln
import pytest


@pytest.fixture
def cellxgene_defaults() -> Generator:
    ln.examples.cellxgene.save_cellxgene_defaults()

    yield

    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    ln.ULabel.filter(type__isnull=False).delete(permanent=True)
    for entity in [
        bt.Disease,
        bt.Ethnicity,
        bt.DevelopmentalStage,
        bt.Phenotype,
        bt.CellType,
        ln.ULabel,
    ]:
        entity.filter().delete(permanent=True)


def test_cellxgene_curation(cellxgene_defaults) -> None:
    """Tests validating a recent CELLxGENE dataset."""
    ln.examples.cellxgene.save_cellxgene_defaults()

    cxg_schema = ln.examples.cellxgene.create_cellxgene_schema(
        field_types="ontology_id",
        organism="mouse",
        spatial_library_id="Thymus_Visium_Exp3A_V2S1_3wk_B6-WT",
    )

    adata = ln.examples.datasets.visium_mouse_cellxgene()

    curator = ln.curators.AnnDataCurator(adata, cxg_schema)
    curator.validate()

    cxg_schema.delete(permanent=True)
