import lamindb as ln
import pandas as pd

union_feature = ln.Feature(
    name="mixed_feature",
    dtype="cat[bionty.Tissue.ontology_id|bionty.CellType.ontology_id]",
).save()

df_mixed = pd.DataFrame({"mixed_feature": ["UBERON:0000178", "CL:0000540"]})

schema = ln.Schema(features=[union_feature], coerce=True).save()

curator = ln.curators.DataFrameCurator(df_mixed, schema)
curator.validate()
