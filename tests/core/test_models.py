import re
import textwrap


def _strip_ansi(text: str) -> str:
    """Remove ANSI escape sequences from a string."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


def test_registry__repr__param():
    import lamindb.models as ln

    param = ln.Param
    expected_repr = textwrap.dedent("""\
    Param
      Simple fields
        .name: CharField
        .dtype: CharField
        .type: CharField
        .created_at: DateTimeField
        .updated_at: DateTimeField
        .aux: JSONField
      Relational fields
        .created_by: User
        .run: Run
        .space: Space
        .values: ParamValue
    """).strip()

    actual_repr = _strip_ansi(repr(param))
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
        .aux: JSONField
      Relational fields
        .space: Space
        .storage: Storage
        .run: Run
        .created_by: User
        .ulabels: ULabel
        .input_of_runs: Run
        .feature_sets: FeatureSet
        .collections: Collection
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
