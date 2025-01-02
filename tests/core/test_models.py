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
        .created_at: DateTimeField
        .updated_at: DateTimeField
      Relational fields
        .created_by: User
        .run: Run
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
        .type: CharField
        .size: BigIntegerField
        .hash: CharField
        .n_objects: BigIntegerField
        .n_observations: BigIntegerField
        .visibility: SmallIntegerField
        .version: CharField
        .is_latest: BooleanField
        .created_at: DateTimeField
        .updated_at: DateTimeField
      Relational fields
        .storage: Storage
        .transform: Transform
        .run: Run
        .created_by: User
        .ulabels: ULabel
        .input_of_runs: Run
        .feature_sets: FeatureSet
        .collections: Collection
    """).strip()

    actual_repr = _strip_ansi(repr(artifact))
    print(actual_repr)
    assert actual_repr.strip() == expected_repr.strip()
