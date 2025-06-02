import lamindb as ln
import bionty as bt
import tiledbsoma as soma
import tiledbsoma.io

adata = ln.core.datasets.small_dataset1(otype="AnnData")
tiledbsoma.io.from_anndata("small_dataset.tiledbsoma", adata, measurement_name="RNA")
experiment = soma.Experiment.open("small_dataset.tiledbsoma")

obs_schema = ln.Schema(
    name="soma_obs_schema",
    features=[
        ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
        ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
    ],
).save()

var_schema = ln.Schema(
    name="soma_var_schema",
    features=[
        ln.Feature(name="var_id", dtype=bt.Gene.ensembl_gene_id).save(),
    ],
    coerce_dtype=True,
).save()

soma_schema = ln.Schema(
    name="soma_experiment_schema",
    otype="tiledbsoma",
    slots={
        "obs": obs_schema,
        "ms:RNA.T": var_schema,
    },
).save()

curator = ln.curators.TiledbsomaExperimentCurator(experiment, soma_schema)
curator.validate()
artifact = curator.save_artifact(
    key="examples/soma_experiment.tiledbsoma",
    description="SOMA experiment with schema validation",
)
assert artifact.schema == soma_schema
artifact.describe()
