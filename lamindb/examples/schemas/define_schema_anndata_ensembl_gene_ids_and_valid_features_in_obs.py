import bionty as bt

import lamindb as ln

cells = ln.Feature(name="Cells", is_type=True).save()
ln.Feature(name="cell_type_by_model", dtype=bt.CellType, type=cells).save()
obs_schema = ln.Schema(itype=cells).save()
var_schema = ln.Schema(
    name="valid_ensembl_gene_ids", itype=bt.Gene.ensembl_gene_id
).save()
schema = ln.Schema(
    name="anndata_ensembl_gene_ids_and_valid_features_in_obs",
    otype="AnnData",
    slots={"obs": obs_schema, "var.T": var_schema},
).save()
