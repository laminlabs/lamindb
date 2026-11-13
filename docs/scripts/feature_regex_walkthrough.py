"""Run with the regex feature branch installed: python feature_regex_walkthrough.py.

Each case defines its own features, schema, and data. There are no assertions.
Expected failures are caught and printed; unexpected exceptions stop the script.
The runner uses a temporary local SQLite database, settings, and storage, then
removes them. It does not connect to playground or your current instance.
"""

import os
from pathlib import Path
from tempfile import TemporaryDirectory


def show_error(error):
    """Print the real exception and any underlying Pandera failure table."""
    print(f"{type(error).__name__}: {error}")
    cause = error
    while cause is not None:
        if hasattr(cause, "failure_cases"):
            print("Underlying Pandera failure cases:")
            print(cause.failure_cases.to_string(index=False))
            break
        cause = cause.__cause__


def valid_identifiers():
    import lamindb as ln
    import pandas as pd

    print("\n=== 1. Valid identifiers ===")
    feature = ln.Feature(name="valid_sample", dtype=str, regex=r"[A-Za-z0-9_-]+").save()
    schema = ln.Schema([feature]).save()
    frame = pd.DataFrame({"valid_sample": ["sample_1", "sample-2", "SAMPLE3"]})
    print("Input:")
    print(frame.to_string())
    print("Expected: validation succeeds; the dataframe is unchanged.")
    print("Actual:")
    ln.curators.DataFrameCurator(frame, schema).validate()
    print("Validation succeeded.")
    print(frame.to_string())


def invalid_identifiers():
    import lamindb as ln
    import pandas as pd
    from lamindb.errors import ValidationError

    print("\n=== 2. Spaces, parentheses, colon, and trailing newline ===")
    feature = ln.Feature(
        name="invalid_sample", dtype=str, regex=r"[A-Za-z0-9_-]+"
    ).save()
    schema = ln.Schema([feature]).save()
    values = ["sample_1", "sample 2", "sample(3)", "sample:4", "sample_5\n"]
    frame = pd.DataFrame({"invalid_sample": values})
    print(f"Input values: {values!r}")
    print("Expected: ValidationError for rows 1, 2, 3, and 4; row 0 passes.")
    print("Actual:")
    try:
        ln.curators.DataFrameCurator(frame, schema).validate()
    except ValidationError as error:
        show_error(error)
    else:
        print("Validation succeeded unexpectedly.")


def wrong_type_and_regex():
    import lamindb as ln
    import pandas as pd
    from lamindb.errors import ValidationError

    print("\n=== 3. Floats supplied to a string feature with a regex ===")
    feature = ln.Feature(name="numeric_sample", dtype=str, regex=r"sample_\d+").save()
    schema = ln.Schema([feature]).save()
    frame = pd.DataFrame({"numeric_sample": [1.5, 2.5]})
    print(frame.to_string())
    print(
        "Expected: one aggregated ValidationError containing both dtype and regex failures."
    )
    print(
        "Expected: numbers are rejected, not converted to strings; no TypeError from re."
    )
    print("Actual:")
    try:
        ln.curators.DataFrameCurator(frame, schema).validate()
    except ValidationError as error:
        show_error(error)
    else:
        print("Validation succeeded unexpectedly.")


def nullability():
    import lamindb as ln
    import pandas as pd
    from lamindb.errors import ValidationError

    print("\n=== 4. Nullability is separate from regex matching ===")
    feature = ln.Feature(
        name="nullable_sample", dtype=str, regex="A+", nullable=True
    ).save()
    schema = ln.Schema([feature]).save()
    frame = pd.DataFrame({"nullable_sample": ["AA", None]})
    print(f"Input values: {frame['nullable_sample'].tolist()!r}")
    print("Expected with nullable=True: validation succeeds; None is skipped by regex.")
    print("Actual:")
    ln.curators.DataFrameCurator(frame, schema).validate()
    print("Validation succeeded.")

    feature.nullable = False
    feature.save()
    print(
        "Expected with nullable=False: nullability failure at row 1, no regex failure."
    )
    print("Actual (using a newly constructed curator):")
    try:
        ln.curators.DataFrameCurator(frame, schema).validate()
    except ValidationError as error:
        show_error(error)
    else:
        print("Validation succeeded unexpectedly.")


def malformed_pattern():
    import lamindb as ln
    from lamindb.errors import ValidationError

    print("\n=== 5. Invalid regex syntax fails before any dataframe validation ===")
    feature = ln.Feature(name="syntax_sample", dtype=str, regex="A+").save()
    print("Expected: assigning '[' raises ValidationError; the old pattern stays 'A+'.")
    print("Actual:")
    try:
        feature.regex = "["
    except ValidationError as error:
        show_error(error)
    else:
        print("Assignment succeeded unexpectedly.")
    print(f"In-memory regex: {feature.regex!r}")
    print(f"Stored regex: {ln.Feature.get(uid=feature.uid).regex!r}")


def constructor_reuse_and_default():
    import lamindb as ln
    from lamindb.errors import ValidationError

    print("\n=== 6. Constructor reuse preserves an omitted default ===")
    original = ln.Feature(
        name="default_sample", dtype=str, regex="A+", default_value="AA"
    ).save()
    print("Initial saved definition: regex='A+', default_value='AA'.")
    print(
        "Expected: reconstructing with regex='B+' and no default raises ValidationError."
    )
    print(
        "Expected: the preserved default 'AA' fails 'B+'; the stored definition is unchanged."
    )
    print("Actual:")
    try:
        ln.Feature(name="default_sample", dtype=str, regex="B+")
    except ValidationError as error:
        show_error(error)
    else:
        print("Construction succeeded unexpectedly.")
    stored = ln.Feature.get(uid=original.uid)
    print(f"Stored regex={stored.regex!r}, default_value={stored.default_value!r}")

    print("Expected: explicitly supplying regex='B+' and default_value='BB' succeeds.")
    print("Expected: constructor changes are in memory until save().")
    print("Actual:")
    replacement = ln.Feature(
        name="default_sample", dtype=str, regex="B+", default_value="BB"
    )
    print(
        f"In memory: regex={replacement.regex!r}, default_value={replacement.default_value!r}"
    )
    stored = ln.Feature.get(uid=original.uid)
    print(
        f"Before save: regex={stored.regex!r}, default_value={stored.default_value!r}"
    )
    replacement.save()
    stored = ln.Feature.get(uid=original.uid)
    print(f"After save: regex={stored.regex!r}, default_value={stored.default_value!r}")


def empty_pattern_and_removal():
    import lamindb as ln
    import pandas as pd
    from lamindb.errors import ValidationError

    print("\n=== 7. An empty pattern is a rule; None removes it ===")
    feature = ln.Feature(name="empty_sample", dtype=str, regex="").save()
    schema = ln.Schema([feature]).save()
    print("Expected: the empty string passes the empty pattern.")
    print("Actual:")
    ln.curators.DataFrameCurator(
        pd.DataFrame({"empty_sample": [""]}), schema
    ).validate()
    print("Validation succeeded.")

    frame = pd.DataFrame({"empty_sample": ["anything"]})
    print("Expected: 'anything' fails the empty pattern.")
    print("Actual:")
    try:
        ln.curators.DataFrameCurator(frame, schema).validate()
    except ValidationError as error:
        show_error(error)
    else:
        print("Validation succeeded unexpectedly.")

    feature.regex = None
    feature.save()
    print("Expected after assigning regex=None and saving: 'anything' passes.")
    print("Actual (using a newly constructed curator):")
    ln.curators.DataFrameCurator(frame, schema).validate()
    print("Validation succeeded.")


def artifact_dataset():
    import lamindb as ln
    import pandas as pd
    from lamindb.errors import ValidationError

    print("\n=== 8. Validate values inside an artifact's dataframe ===")
    feature = ln.Feature(name="dataset_sample", dtype=str, regex=r"sample_\d+").save()
    schema = ln.Schema([feature]).save()
    print(
        "Expected: Artifact.from_dataframe(..., schema=...) rejects 'sample(bad)' before save()."
    )
    print("Actual:")
    try:
        ln.Artifact.from_dataframe(
            pd.DataFrame({"dataset_sample": ["sample(bad)"]}),
            key="regex-demo/invalid-dataset.parquet",
            schema=schema,
        )
    except ValidationError as error:
        show_error(error)
    else:
        print("Artifact construction succeeded unexpectedly.")

    print(
        "Expected: a dataframe containing 'sample_1' validates and saves successfully."
    )
    print("Actual:")
    artifact = ln.Artifact.from_dataframe(
        pd.DataFrame({"dataset_sample": ["sample_1"]}),
        key="regex-demo/valid-dataset.parquet",
        schema=schema,
    ).save()
    print(f"Saved artifact: {artifact.key}")
    print(artifact.load().to_string())


def artifact_annotation():
    import lamindb as ln
    import pandas as pd
    from lamindb.errors import ValidationError

    print("\n=== 9. Validate an artifact annotation, not its dataframe columns ===")
    feature = ln.Feature(
        name="annotation_sample", dtype=str, regex=r"sample_\d+"
    ).save()
    artifact = ln.Artifact.from_dataframe(
        pd.DataFrame({"measurement": [10.0]}),
        key="regex-demo/annotated-dataset.parquet",
    ).save()
    print("Expected: setting annotation_sample='sample_1' succeeds.")
    print("Actual:")
    artifact.features.set_values({feature: "sample_1"})
    print(f"Stored annotation: {artifact.features.get_values()[feature.name]!r}")

    print("Expected: updating the annotation to 'sample(bad)' raises ValidationError.")
    print(
        "Expected: the error includes artifact/assignment context; stored value stays 'sample_1'."
    )
    print("Actual:")
    try:
        artifact.features.update_values({feature: "sample(bad)"})
    except ValidationError as error:
        show_error(error)
    else:
        print("Annotation update succeeded unexpectedly.")
    print(
        f"Stored annotation after rejected update: {artifact.features.get_values()[feature.name]!r}"
    )


def main():
    """Infrastructure only: each case above owns its feature/schema/data setup."""
    original_directory = Path.cwd()
    with TemporaryDirectory(prefix="lamindb-regex-demo-") as directory:
        # Set this before importing LaminDB so existing account settings are untouched.
        os.environ["LAMIN_SETTINGS_DIR"] = directory
        os.environ["LAMIN_CACHE_DIR"] = str(Path(directory) / "cache")
        os.chdir(directory)
        try:
            import lamindb as ln
            from django.db import connections
            from lamindb_setup import logger

            # Keep the output focused on the explicit expected/actual prints.
            with logger.mute():
                ln.setup.init(storage=Path(directory) / "storage", name="regex-demo")
                try:
                    print(
                        "Regex walkthrough: temporary local SQLite database; no assertions."
                    )
                    valid_identifiers()
                    invalid_identifiers()
                    wrong_type_and_regex()
                    nullability()
                    malformed_pattern()
                    constructor_reuse_and_default()
                    empty_pattern_and_removal()
                    artifact_dataset()
                    artifact_annotation()
                    print(
                        "\nAll scenarios finished. Temporary data will now be removed."
                    )
                finally:
                    connections.close_all()
        finally:
            os.chdir(original_directory)


if __name__ == "__main__":
    main()
