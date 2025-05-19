import lamindb as ln
import bionty as bt


# define the global obs schema
obs_schema = ln.Schema(
    name="mudata_papalexi21_subset_obs_schema",
    features=[
        ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
        ln.Feature(name="replicate", dtype="cat[ULabel[Replicate]]").save(),
    ],
).save()

# define the ['rna'].obs schema
obs_schema_rna = ln.Schema(
    name="mudata_papalexi21_subset_rna_obs_schema",
    features=[
        ln.Feature(name="nCount_RNA", dtype=int).save(),
        ln.Feature(name="nFeature_RNA", dtype=int).save(),
        ln.Feature(name="percent.mito", dtype=float).save(),
    ],
).save()

# define the ['hto'].obs schema
obs_schema_hto = ln.Schema(
    name="mudata_papalexi21_subset_hto_obs_schema",
    features=[
        ln.Feature(name="nCount_HTO", dtype=int).save(),
        ln.Feature(name="nFeature_HTO", dtype=int).save(),
        ln.Feature(name="technique", dtype=bt.ExperimentalFactor).save(),
    ],
).save()

# define ['rna'].var schema
var_schema_rna = ln.Schema(
    name="mudata_papalexi21_subset_rna_var_schema",
    itype=bt.Gene.symbol,
    dtype=float,
).save()

# define composite schema
mudata_schema = ln.Schema(
    name="mudata_papalexi21_subset_mudata_schema",
    otype="MuData",
    slots={
        "obs": obs_schema,
        "rna:obs": obs_schema_rna,
        "hto:obs": obs_schema_hto,
        "rna:var": var_schema_rna,
    },
).save()

# curate a MuData
mdata = ln.core.datasets.mudata_papalexi21_subset()
bt.settings.organism = "human"  # set the organism to map gene symbols
curator = ln.curators.MuDataCurator(mdata, mudata_schema)
artifact = curator.save_artifact(key="examples/mudata_papalexi21_subset.h5mu")
assert artifact.schema == mudata_schema
