from pathlib import Path

import lamindb as ln
from _dataset_fixtures import (  # noqa
    get_mini_csv,
)

# Without escaping, rich parses `[link=...]` into an unescaped `href` in Jupyter
# HTML and raises `MarkupError` on the unmatched `[/b]` closing tag.
PAYLOAD = '[link="><b>]x[/b]'


def test_describe_shows_db_strings_literally(get_mini_csv: Path):  # noqa: F811
    label = ln.Record(name=PAYLOAD).save()
    feature = ln.Feature(name="describe_markup_value", dtype=str).save()
    artifact = ln.Artifact(get_mini_csv, description="describe markup").save()
    artifact.records.add(label)
    artifact.features.add_values({feature.name: PAYLOAD})
    member = ln.Feature(name="dmm", dtype=str, default_value=PAYLOAD).save()
    schema = ln.Schema([member], name="describe_markup_schema").save()
    transform = ln.Transform(key="describe_markup.py").save()
    run = ln.Run(transform, params={PAYLOAD: PAYLOAD}).save()
    try:
        # once in Entities (labels), once in Features (external values)
        assert artifact.describe(return_str=True).count(PAYLOAD) == 2
        assert PAYLOAD in schema.describe(return_str=True)
        assert f"{PAYLOAD}: {PAYLOAD}" in run.describe(return_str=True)
    finally:
        run.delete(permanent=True)
        transform.delete(permanent=True)
        artifact.delete(permanent=True)
        schema.delete(permanent=True)
        member.delete(permanent=True)
        feature.delete(permanent=True)
        label.delete(permanent=True)
