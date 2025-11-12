import textwrap

import lamindb as ln
import pandas as pd
import pytest
from lamindb.models._describe import strip_ansi_from_string as _strip_ansi


def test_registry__repr__feature():
    import lamindb.models as ln

    feature = ln.Feature
    expected_repr = textwrap.dedent("""\
    Feature
      Simple fields
        .uid: CharField
        .name: CharField
        .dtype: CharField
        .is_type: BooleanField
        .unit: CharField
        .description: TextField
        .array_rank: SmallIntegerField
        .array_size: IntegerField
        .array_shape: JSONField
        .proxy_dtype: CharField
        .synonyms: TextField
        .is_locked: BooleanField
        .created_at: DateTimeField
        .updated_at: DateTimeField
      Relational fields
        .branch: Branch
        .space: Space
        .created_by: User
        .run: Run
        .type: Feature
        .schemas: Schema
        .features: Feature
        .values: FeatureValue
        .projects: Project
        .blocks: FeatureBlock
    """).strip()

    actual_repr = _strip_ansi(repr(feature))
    print(actual_repr)
    assert actual_repr.strip() == expected_repr.strip()


def test_registry__repr__artifact():
    import lamindb.models as ln

    artifact = ln.Artifact
    expected_repr = textwrap.dedent("""\
    Artifact
      Simple fields
        .uid: CharField
        .key: CharField
        .description: TextField
        .suffix: CharField
        .kind: CharField
        .otype: CharField
        .size: BigIntegerField
        .hash: CharField
        .n_files: BigIntegerField
        .n_observations: BigIntegerField
        .version: CharField
        .is_latest: BooleanField
        .is_locked: BooleanField
        .created_at: DateTimeField
        .updated_at: DateTimeField
      Relational fields
        .branch: Branch
        .space: Space
        .storage: Storage
        .run: Run
        .schema: Schema
        .created_by: User
        .input_of_runs: Run
        .feature_sets: Schema
        .linked_in_records: Record
        .users: User
        .ulabels: ULabel
        .collections: Collection
        .records: Record
        .references: Reference
        .projects: Project
        .blocks: Block
      Bionty fields
        .organisms: bionty.Organism
        .genes: bionty.Gene
        .proteins: bionty.Protein
        .cell_markers: bionty.CellMarker
        .tissues: bionty.Tissue
        .cell_types: bionty.CellType
        .diseases: bionty.Disease
        .cell_lines: bionty.CellLine
        .phenotypes: bionty.Phenotype
        .pathways: bionty.Pathway
        .experimental_factors: bionty.ExperimentalFactor
        .developmental_stages: bionty.DevelopmentalStage
        .ethnicities: bionty.Ethnicity
    """).strip()

    actual_repr = _strip_ansi(repr(artifact))
    print(actual_repr)
    assert actual_repr.strip() == expected_repr.strip()


def test_unsaved_relationship_modification_attempts():
    af = ln.Artifact.from_dataframe(
        pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}), description="testme"
    )

    new_label = ln.Record(name="testlabel").save()
    with pytest.raises(ValueError) as excinfo:
        af.records.add(new_label)

    assert (
        str(excinfo.value)
        == "You are trying to access the many-to-many relationships of an unsaved Artifact object. Please save it first using '.save()'."
    )

    new_label.delete(permanent=True)
    af.delete(permanent=True)


def test_unsaved_model_different_instance():
    af = ln.Artifact.connect("laminlabs/lamindata").get(
        key="scrna/micro-macfarland2020.h5ad"
    )

    new_label = ln.Record(name="testlabel").save()
    with pytest.raises(ValueError) as excinfo:
        af.records.add(new_label)

    assert (
        str(excinfo.value)
        == "Cannot label a record from instance 'laminlabs/lamindata'. "
        "Please save the record first to your instance using '.save()'."
    )

    new_label.delete(permanent=True)
