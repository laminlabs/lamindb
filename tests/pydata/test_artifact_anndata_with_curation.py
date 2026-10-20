import lamindb as ln


def test_create_anndata_with_curation():
    adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
    schema = ln.examples.schemas.anndata_ensembl_gene_ids_and_obs_feature_type()
    cells = ln.Feature.get(name="Cells", is_type=True)
    feature1 = ln.Feature(name="sample_note", dtype=str, type=cells).save()

    # ingest the first time
    artifact = ln.Artifact.from_anndata(
        adata,
        key="examples/mini_immuno1.h5ad",
        schema=schema,
    ).save()
    # capture the obs_schema because we'll overwrite it
    obs_schema = artifact.features.slots["obs"]

    # define another feature so that upon re-ingestion, we track more than before
    # (this also tests non-trivial idempotency)
    feature2 = ln.Feature(name="treatment_time_h", dtype=int, type=cells).save()
    artifact = ln.Artifact.from_anndata(
        adata,
        key="examples/mini_immuno1.h5ad",
        schema=schema,
    ).save()

    schemas = artifact.features.slots
    artifact.delete(permanent=True)
    for schema in schemas.values():
        schema.delete(permanent=True)
    obs_schema.delete(permanent=True)
    feature1.delete(permanent=True)
    feature2.delete(permanent=True)
