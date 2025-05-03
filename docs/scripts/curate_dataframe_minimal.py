import lamindb as ln
import bionty as bt

# define valid labels
perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
ln.ULabel(name="DMSO", type=perturbation).save()
ln.ULabel(name="IFNG", type=perturbation).save()
bt.CellType.from_source(name="B cell").save()
bt.CellType.from_source(name="T cell").save()

# define schema
schema = ln.Schema(
    name="small_dataset1_obs_level_metadata",
    features=[
        ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
        ln.Feature(name="sample_note", dtype=str).save(),
        ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
        ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
    ],
).save()

# curate a DataFrame
df = ln.core.datasets.small_dataset1(otype="DataFrame")
curator = ln.curators.DataFrameCurator(df, schema)
artifact = curator.save_artifact(key="examples/dataset1.parquet")
assert artifact.schema == schema
