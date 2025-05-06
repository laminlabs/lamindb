import lamindb as ln
import bionty as bt

obs_schema = ln.examples.schemas.valid_features()
varT_schema = ln.Schema(
    name="valid_ensembl_gene_ids", itype=bt.Gene.ensembl_gene_id
).save()
schema = ln.Schema(
    name="anndata_ensembl_gene_ids_and_valid_features_in_obs",
    otype="AnnData",
    slots={"obs": obs_schema, "var.T": varT_schema},
).save()
