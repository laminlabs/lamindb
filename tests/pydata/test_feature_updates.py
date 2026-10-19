import lamindb as ln
import pytest


@pytest.fixture(params=["Record", "Artifact", "Run"])
def annotations(request):
    score = ln.Feature(name="update_score", dtype=int).save()
    numbers = ln.Feature(name="update_numbers", dtype=list[int]).save()
    note = ln.Feature(name="update_note", dtype=str).save()
    category = ln.Feature(name="update_category", dtype=ln.Record).save()
    categories = ln.Feature(name="update_categories", dtype=list[ln.Record]).save()
    labels = [ln.Record(name=f"update_label_{i}").save() for i in range(3)]
    transform = None
    if request.param == "Record":
        host = ln.Record(name="update_record").save()
    elif request.param == "Artifact":
        host = ln.Artifact(
            ".gitignore", key="update_artifact", skip_hash_lookup=True
        ).save()
    else:
        transform = ln.Transform(key="update_transform").save()
        host = ln.Run(transform, name="update_run").save()
    features = [score, numbers, note, category, categories]
    try:
        yield host, features, labels
    finally:
        host.delete(permanent=True)
        if transform is not None:
            transform.delete(permanent=True)
        for label in labels:
            label.delete(permanent=True)
        for feature in features:
            feature.delete(permanent=True)


def test_update_values_replaces_only_supplied_features(annotations):
    host, (score, numbers, note, category, categories), labels = annotations
    host.features.set_values(
        {
            score: 1,
            numbers: [1, 2],
            note: "keep",
            category: labels[0],
            categories: labels[:2],
        }
    )
    # Untouched annotation rows retain their identity as well as their value.
    if isinstance(host, ln.Record):
        note_id = host.values_json.get(feature=note).id
    else:
        note_id = host.json_values.through.objects.get(
            **{host.__class__.__name__.lower(): host, "jsonvalue__feature": note}
        ).id

    host.features.update_values(
        {score: 2, numbers.name: [3], category: labels[1].name, categories: labels[1:]}
    )
    assert host.features.get_values() == {
        score.name: 2,
        numbers.name: [3],
        note.name: "keep",
        category.name: labels[1].name,
        categories.name: [label.name for label in labels[1:]],
    }
    if isinstance(host, ln.Record):
        assert host.values_json.get(feature=note).id == note_id
    else:
        assert (
            host.json_values.through.objects.get(
                **{host.__class__.__name__.lower(): host, "jsonvalue__feature": note}
            ).id
            == note_id
        )

    host.features.update_values({})
    host.features.update_values({numbers: [], category: None, categories: None})
    assert host.features.get_values() == {
        score.name: 2,
        numbers.name: [],
        note.name: "keep",
    }
    host.features.update_values({score: None, numbers: None})
    assert host.features.get_values() == {note.name: "keep"}


def test_update_values_adds_missing_and_validates_before_removing(annotations):
    host, (score, numbers, note, category, categories), labels = annotations
    host.features.set_values({score: 1, categories: labels[:2]})
    # A missing feature first must not stop replacement of later features.
    host.features.update_values(
        {note: "new", category: labels[2], score: 2, categories: labels[1:]}
    )
    before = host.features.get_values()
    with pytest.raises(ln.errors.ValidationError):
        host.features.update_values({score: "invalid", categories: [labels[0]]})
    assert host.features.get_values() == before
    with pytest.raises(ln.errors.ValidationError, match="Conflicting values"):
        host.features.update_values({score: 3, score.name: 4})
    unsaved = ln.Record(name="update_unsaved")
    with pytest.raises(ln.errors.ValidationError, match="not saved"):
        host.features.update_values({category: unsaved})
    assert host.features.get_values() == before


def test_update_values_rolls_back_failed_insertion(annotations, monkeypatch):
    from lamindb.models import FeatureManager

    host, (score, numbers, note, category, categories), labels = annotations
    host.features.set_values({score: 1, categories: labels[:2]})
    before = host.features.get_values()
    original = FeatureManager._add_values

    def fail_after_insert(self, *args, **kwargs):
        original(self, *args, **kwargs)
        raise RuntimeError("insertion failed")

    monkeypatch.setattr(FeatureManager, "_add_values", fail_after_insert)
    with pytest.raises(RuntimeError, match="insertion failed"):
        host.features.update_values({score: 2, categories: labels[1:]})
    assert host.features.get_values() == before


@pytest.mark.parametrize("inherited", [False, True])
def test_update_values_validates_complete_schema(inherited):
    score = ln.Feature(name="update_schema_score", dtype=int, nullable=False).save()
    note = ln.Feature(name="update_schema_note", dtype=str).save()
    index = ln.Feature(name="update_schema_index", dtype=str).save()
    schema = ln.Schema([score, note.with_config(optional=True)], index=index).save()
    frame = ln.Record(name="update_schema_frame", is_type=True, schema=schema).save()
    kwargs = {"type": frame} if inherited else {}
    record = ln.Record(
        name="row", features={score: 1, note: "keep", index: "row"}, **kwargs
    ).save()
    schema_kwargs = {} if inherited else {"schema": schema}
    try:
        record.features.update_values({score: 2}, **schema_kwargs)
        assert record.features.get_values()[note.name] == "keep"
        assert record.name == "row"
        record.features.update_values({note: None}, **schema_kwargs)
        assert note.name not in record.features.get_values()
        with pytest.raises(ln.errors.ValidationError):
            record.features.update_values({score: None}, **schema_kwargs)
        assert record.features.get_values()[score.name] == 2
        record.features.update_values({index: "renamed"}, **schema_kwargs)
        assert ln.Record.get(id=record.id).name == "renamed"
        assert record.features.get_values()[score.name] == 2
        assert record.features.update_values({}) is None
    finally:
        record.delete(permanent=True)
        frame.delete(permanent=True)
        schema.delete(permanent=True)
        for feature in (score, note, index):
            feature.delete(permanent=True)


def test_update_values_preserves_ambiguous_categorical_identity():
    value_type = ln.Record(name="update_value_type", is_type=True).save()
    a = ln.Record(name="update_shared_name", type=value_type).save()
    b = ln.Record(name="update_shared_name", type=value_type).save()
    category = ln.Feature(name="update_identity", dtype=value_type).save()
    score = ln.Feature(name="update_identity_score", dtype=int).save()
    schema = ln.Schema([category, score]).save()
    record = ln.Record(name="update_identity_record").save()
    try:
        record.features.set_values({category: a, score: 1}, schema=schema)
        link = record.values_record.get(feature=category)
        record.features.update_values({score: 2}, schema=schema)
        assert record.values_record.get(feature=category).id == link.id
        assert record.values_record.get(feature=category).value_id == a.id
        assert record.features.get_values()[score.name] == 2
    finally:
        record.delete(permanent=True)
        a.delete(permanent=True)
        b.delete(permanent=True)
        schema.delete(permanent=True)
        category.delete(permanent=True)
        score.delete(permanent=True)
        value_type.delete(permanent=True)


def test_update_values_mapped_fields_and_inverse_relationships():
    author = ln.Feature(name="update_author", dtype=ln.Record).save()
    books = ln.Feature(
        name="update_books",
        dtype=list[ln.Record],
        values_through=author,
    ).save()
    description = ln.Feature(
        name="update_description", dtype=str, values_through="description"
    ).save()
    score = ln.Feature(name="update_author_score", dtype=int).save()
    alice = ln.Record(name="update_alice").save()
    bob = ln.Record(name="update_bob").save()
    first = ln.Record(name="update_first_book").save()
    second = ln.Record(name="update_second_book").save()
    try:
        alice.features.set_values({books: [first], description: "old", score: 1})
        bob.features.set_values({books: [second]})
        alice.features.update_values({books: [second], description: "new"})
        assert not first.values_record.filter(feature=author).exists()
        assert second.values_record.get(feature=author).value_id == alice.id
        assert alice.features.get_values()[score.name] == 1
        assert ln.Record.get(id=alice.id).description == "new"
        alice.features.update_values({description: None, books: None})
        assert not second.values_record.filter(feature=author).exists()
        assert ln.Record.get(id=alice.id).description is None
        assert alice.features.get_values()[score.name] == 1
        assert alice.features.get_values().get(books.name) in (None, [])
    finally:
        for record in (first, second, alice, bob):
            record.delete(permanent=True)
        for feature in (books, author, description, score):
            feature.delete(permanent=True)


def test_update_values_complete_schema_across_hosts(annotations):
    host, (score, numbers, note, category, categories), labels = annotations
    schema = ln.Schema([score, category, note.with_config(optional=True)]).save()
    kwargs = {"schema": schema}
    if isinstance(host, ln.Artifact):
        host.schema = schema
        host.save()
        kwargs = {}
    try:
        host.features.set_values(
            {score: 1, category: labels[0], note: "keep"}, **kwargs
        )
        category_link_id = host.features[category.name].id
        host.features.update_values({score: 2}, **kwargs)
        assert host.features.get_values()[note.name] == "keep"
        assert host.features[category.name].id == category_link_id
        host.features.update_values({note: None}, **kwargs)
        assert note.name not in host.features.get_values()
        with pytest.raises(ln.errors.ValidationError):
            host.features.update_values({category: None}, **kwargs)
        assert host.features[category.name] == labels[0]
        assert host.features.get_values()[score.name] == 2
    finally:
        if isinstance(host, ln.Artifact):
            host.schema = None
            host.save()
        schema.delete(permanent=True)


def test_update_values_protects_artifact_dataset_features():
    score = ln.Feature(name="update_dataset_score", dtype=int).save()
    note = ln.Feature(name="update_external_note", dtype=str).save()
    schema = ln.Schema([score]).save()
    artifact = ln.Artifact(
        ".gitignore", key="update_dataset_protection", skip_hash_lookup=True
    ).save()
    try:
        artifact.features.set_values({score: 1, note: "keep"})
        artifact.features._add_schema(schema, "columns")
        with pytest.raises(
            ValueError, match="Cannot remove values for dataset features"
        ):
            artifact.features.update_values({note: "changed", score: 2})
        assert artifact.features.get_values() == {score.name: 1, note.name: "keep"}
    finally:
        artifact.delete(permanent=True)
        schema.delete(permanent=True)
        score.delete(permanent=True)
        note.delete(permanent=True)


def test_update_values_with_duplicate_feature_names():
    feature_type = ln.Feature(name="update_feature_scope", is_type=True).save()
    global_feature = ln.Feature(name="update_scoped_value", dtype=str).save()
    typed_feature = ln.Feature(
        name="update_scoped_value", dtype=str, type=feature_type
    ).save()
    record = ln.Record(name="update_scoped_record").save()
    try:
        record.features.set_values({global_feature: "global"})
        record.features.update_values({typed_feature: "typed"})
        assert record.values_json.get(feature=global_feature).value == "global"
        assert record.values_json.get(feature=typed_feature).value == "typed"
        record.features.update_values({typed_feature: "new"})
        assert record.values_json.get(feature=global_feature).value == "global"
        assert record.values_json.get(feature=typed_feature).value == "new"
    finally:
        record.delete(permanent=True)
        typed_feature.delete(permanent=True)
        global_feature.delete(permanent=True)
        feature_type.delete(permanent=True)


def test_update_values_handles_duplicate_categorical_links(annotations):
    host, (score, numbers, note, category, categories), labels = annotations
    host.features.set_values({score: 1, categories: labels[:2]})
    host.features.update_values({categories: [labels[2], labels[0], labels[2]]})
    assert host.features.get_values()[categories.name] == [
        labels[2].name,
        labels[0].name,
    ]
    assert host.features.get_values()[score.name] == 1


def test_update_values_replaces_all_registries_of_categorical_union():
    feature = ln.Feature(name="update_union", dtype="list[cat[Record|ULabel]]").save()
    record = ln.Record(name="update_union_host").save()
    label = ln.ULabel(name="update_union_label").save()
    value = ln.Record(name="update_union_value").save()
    try:
        record.features.set_values({feature: [value, label]})
        record.features.update_values({feature: [value]})
        assert record.features.get_values()[feature.name] == [value.name]
        assert not record.values_ulabel.filter(feature=feature).exists()
        record.features.update_values({feature: None})
        assert feature.name not in record.features.get_values()
    finally:
        record.delete(permanent=True)
        value.delete(permanent=True)
        label.delete(permanent=True)
        feature.delete(permanent=True)
