# %% [markdown]
# # Curate `AnnData` based on the CELLxGENE schema
#
# This guide shows how to curate an AnnData object with the help of [`laminlabs/cellxgene`](https://lamin.ai/laminlabs/cellxgene) against the [CELLxGENE schema v5.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md).

# %% [markdown]
# Load your instance where you want to register the curated AnnData object:

# %% tags=["hide-output"]
# pip install 'lamindb[bionty,jupyter]' cellxgene-lamin
# cellxgene-schema has pinned dependencies. Therefore we recommend installing it into a separate environment using `uv` or `pipx`
# uv tool install cellxgene-schema==5.1.0
# !lamin init --storage ./test-cellxgene-curate --modules bionty

# %%
import lamindb as ln


def get_semi_curated_dataset():
    adata = ln.core.datasets.anndata_human_immune_cells()
    adata.obs["sex_ontology_term_id"] = "PATO:0000384"
    adata.obs["organism"] = "human"
    adata.obs["sex"] = "unknown"
    # create some typos in the metadata
    adata.obs["tissue"] = adata.obs["tissue"].cat.rename_categories({"lung": "lungg"})
    # new donor ids
    adata.obs["donor"] = adata.obs["donor"].astype(str) + "-1"
    # drop animal cell
    adata = adata[adata.obs["cell_type"] != "animal cell", :]
    # remove columns that are reserved in the cellxgene schema
    adata.var.drop(columns=["feature_reference", "feature_biotype"], inplace=True)
    adata.raw.var.drop(
        columns=["feature_name", "feature_reference", "feature_biotype"], inplace=True
    )
    print("loaded dataset")
    return adata


def test_cxg_curator():
    # %% [markdown]
    # Let's start with an AnnData object that we'd like to inspect and curate.
    # We are writing it to disk to run [CZI's cellxgene-schema CLI tool](https://github.com/chanzuckerberg/single-cell-curation) which verifies whether an on-disk h5ad dataset adheres to the cellxgene schema.

    # %% tags=["hide-output"]
    adata = get_semi_curated_dataset()

    # %% [markdown]
    # Initially, the cellxgene-schema validator of CZI does not pass and we need to curate the dataset.
    # cellxgene-schema validate anndata_human_immune_cells.h5ad

    # %% [markdown]
    # ## Validate and curate metadata

    # %% tags=["hide-output"]
    curator = ln.curators.CellxGeneAnnDataCurator(
        adata, organism="human", schema_version="5.1.0"
    )

    # %% [markdown]
    # Let's fix the "donor_id" column name:

    # %%
    adata.obs.rename(columns={"donor": "donor_id"}, inplace=True)

    # %%
    curator.validate()

    # %% [markdown]
    # For the missing columns, we can pass default values suggested from CELLxGENE which will automatically add them to the AnnData object:

    # %% [markdown]
    # ```{note}
    # CELLxGENE requires columns `tissue`, `organism`, and `assay` to have existing values from the ontologies.
    # Therefore, these columns need to be added and populated manually.
    # ```

    # %% tags=["hide-output"]
    curator = ln.curators.CellxGeneAnnDataCurator(
        adata,
        defaults=ln.curators.CellxGeneFields.OBS_FIELD_DEFAULTS,
        organism="human",
        schema_version="5.1.0",
    )

    # %% tags=["hide-output"]
    assert curator.validate() is not True

    # %% [markdown]
    # ## Remove unvalidated values

    # %% [markdown]
    # We remove all unvalidated genes.
    # These genes may exist in a different release of ensembl but are not valid for the ensembl version of cellxgene schema 5.0.0 (ensembl release 110).

    # %%
    # adata = adata[:, ~adata.var.index.isin(curator.non_validated["var_index"])].copy()
    # if adata.raw is not None:
    #     raw_data = adata.raw.to_adata()
    #     raw_data = raw_data[
    #         :, ~raw_data.var_names.isin(curator.non_validated["var_index"])
    #     ].copy()
    #     adata.raw = raw_data

    # %%
    curator = ln.curators.CellxGeneAnnDataCurator(
        adata, organism="human", schema_version="5.1.0"
    )

    # %% [markdown]
    # ## Register new metadata labels
    #
    # Following the suggestions above to register genes and labels that aren't present in the current instance:
    #
    # (Note that our instance is rather empty. Once you filled up the registries, registering new labels won't be frequently needed)

    # %% [markdown]
    # For donors, we register the new labels:

    # %% tags=["hide-output"]
    curator.add_new_from("donor_id")

    # %% [markdown]
    # An error is shown for the tissue label "lungg", which is a typo, should be "lung". Let's fix it:

    # %% tags=["hide-output"]
    tissues = curator.lookup().tissue

    # %% tags=["hide-output"]
    adata.obs["tissue"] = adata.obs["tissue"].cat.rename_categories(
        {"lungg": tissues.lung.name}
    )

    # %% [markdown]
    # Let's validate the object again:

    # %% tags=["hide-output"]
    assert curator.validate() is True

    # %% tags=["hide-output"]
    adata.obs.head()

    # %% [markdown]
    # ## Save artifact

    # %% tags=["hide-output"]
    artifact = curator.save_artifact(
        description=f"dataset curated against cellxgene schema {curator.schema_version}"
    )

    # %% tags=["hide-output"]
    artifact.describe()

    # %% [markdown]
    # ## Return an input h5ad file for cellxgene-schema

    # %% tags=["hide-output"]
    title = "Cross-tissue immune cell analysis reveals tissue-specific features in humans (for test demo only)"
    curator.to_cellxgene_anndata(is_primary_data=True, title=title)

    # %% [markdown]
    # Now it passes:
    # cellxgene-schema validate anndata_human_immune_cells_cxg.h5ad
    # ```{note}
    #
    # The Curate class is designed to validate all metadata for adherence to ontologies.
    # It does not reimplement all rules of the cellxgene schema and we therefore recommend running the [cellxgene-schema](https://github.com/chanzuckerberg/single-cell-curation) if full adherence beyond metadata is a necessity.
    # ```


test_cxg_curator()
