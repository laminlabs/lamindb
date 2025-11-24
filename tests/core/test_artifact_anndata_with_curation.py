import lamindb as ln


def test_create_anndata_with_curation():
    adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
    feature1 = ln.Feature(name="sample_note", dtype=str).save()

    # ingest the first time
    artifact = ln.Artifact.from_anndata(
        adata,
        key="examples/mini_immuno1.h5ad",
        schema="ensembl_gene_ids_and_valid_features_in_obs",
    ).save()

    # define another feature so that upon re-ingestion, we track more than before
    # (this also tests non-trivial idempotency)
    feature2 = ln.Feature(name="treatment_time_h", dtype=int).save()
    artifact = ln.Artifact.from_anndata(
        adata,
        key="examples/mini_immuno1.h5ad",
        schema="ensembl_gene_ids_and_valid_features_in_obs",
    ).save()

    schemas = artifact.feature_sets.all()
    artifact.delete(permanent=True)
    schemas.delete(permanent=True)
    feature1.delete(permanent=True)
    feature2.delete(permanent=True)
