import bionty as bt

import lamindb as ln

# define valid labels
perturbations = ln.Record(name="Perturbations", is_type=True).save()
ln.Record(name="DMSO", type=perturbations).save()
ln.Record(name="IFNG", type=perturbations).save()
bt.CellType.from_source(name="B cell").save()
bt.CellType.from_source(name="T cell").save()

# define valid features for immuno samples
immuno_features = ln.Feature(name="ImmunoSamples", is_type=True).save()
ln.Feature(name="perturbation", dtype=perturbations, type=immuno_features).save()
ln.Feature(name="cell_type_by_expert", dtype=bt.CellType, type=immuno_features).save()
ln.Feature(name="cell_type_by_model", dtype=bt.CellType, type=immuno_features).save()
ln.Feature(
    name="assay_oid", dtype=bt.ExperimentalFactor.ontology_id, type=immuno_features
).save()
ln.Feature(name="concentration", dtype=str, type=immuno_features).save()
ln.Feature(
    name="treatment_time_h", dtype="num", coerce=True, type=immuno_features
).save()
ln.Feature(name="donor", dtype=str, nullable=True, type=immuno_features).save()
ln.Feature(
    name="donor_ethnicity", dtype=list[bt.Ethnicity], type=immuno_features
).save()
