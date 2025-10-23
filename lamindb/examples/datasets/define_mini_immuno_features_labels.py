import bionty as bt

import lamindb as ln

# define valid labels
perturbation_type = ln.Record(name="Perturbation", is_type=True).save()
ln.Record(name="DMSO", type=perturbation_type).save()
ln.Record(name="IFNG", type=perturbation_type).save()
bt.CellType.from_source(name="B cell").save()
bt.CellType.from_source(name="T cell").save()

# define valid features
ln.Feature(name="perturbation", dtype=perturbation_type).save()
ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save()
ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save()
ln.Feature(name="assay_oid", dtype=bt.ExperimentalFactor.ontology_id).save()
ln.Feature(name="concentration", dtype=str).save()
ln.Feature(name="treatment_time_h", dtype="num", coerce_dtype=True).save()
ln.Feature(name="donor", dtype=str, nullable=True).save()
ln.Feature(name="donor_ethnicity", dtype=list[bt.Ethnicity]).save()
