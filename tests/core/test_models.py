import re
import textwrap

import lamindb as ln
import pandas as pd
import pytest


def _strip_ansi(text: str) -> str:
    """Remove ANSI escape sequences from a string."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


def test_registry__repr__feature():
    import lamindb.models as ln

    feature = ln.Param
    expected_repr = textwrap.dedent("""\
    Feature
      Simple fields
        .uid: CharField
        .name: CharField
        .dtype: CharField
        .is_type: BooleanField
        .unit: CharField
        .description: CharField
        .array_rank: SmallIntegerField
        .array_size: IntegerField
        .array_shape: JSONField
        .proxy_dtype: CharField
        .synonyms: TextField
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
        .description: CharField
        .suffix: CharField
        .kind: CharField
        .otype: CharField
        .size: BigIntegerField
        .hash: CharField
        .n_files: BigIntegerField
        .n_observations: BigIntegerField
        .version: CharField
        .is_latest: BooleanField
        .created_at: DateTimeField
        .updated_at: DateTimeField
      Relational fields
        .branch: Branch
        .space: Space
        .storage: Storage
        .run: Run
        .schema: Schema
        .created_by: User
        .ulabels: ULabel
        .input_of_runs: Run
        .feature_sets: Schema
        .collections: Collection
        .records: Record
        .references: Reference
        .projects: Project
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
    af = ln.Artifact.from_df(
        pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}), description="testme"
    )

    new_label = ln.ULabel(name="testlabel").save()
    with pytest.raises(ValueError) as excinfo:
        af.ulabels.add(new_label)

    assert (
        str(excinfo.value)
        == "You are trying to access the many-to-many relationships of an unsaved Artifact object. Please save it first using '.save()'."
    )

    new_label.delete()
    af.delete()
