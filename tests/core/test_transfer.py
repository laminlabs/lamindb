import lamindb as ln


def test_schema_transfer_defaults_to_annotations():
    db = ln.DB("laminlabs/lamindata")
    schema_uid = "pnQvQVcQ417bfmVq"
    remote_schema = db.Schema.get(schema_uid)
    remote_member_names = remote_schema.members.to_list("name")
    assert len(remote_member_names) > 0
    remote_perturbation_type_name = "Perturbation"

    existing_local = ln.Schema.filter(uid=schema_uid).one_or_none()
    if existing_local is not None:
        existing_local.delete(permanent=True)

    existing_local_perturbation_type = ln.Record.filter(
        name=remote_perturbation_type_name, is_type=True
    ).one_or_none()
    if existing_local_perturbation_type is not None:
        ln.Record.filter(type=existing_local_perturbation_type).delete(permanent=True)
        existing_local_perturbation_type.delete(permanent=True)

    assert (
        ln.Record.filter(name=remote_perturbation_type_name, is_type=True).count() == 0
    )
    assert ln.Record.filter(type__name=remote_perturbation_type_name).count() == 0

    record_only = db.Schema.get(schema_uid).save(transfer="record")
    assert record_only.members.count() == 0

    transferred = db.Schema.get(schema_uid).save()
    transferred_member_names = transferred.members.to_list("name")
    assert transferred_member_names == remote_member_names

    perturbation_feature = transferred.members.get(name="perturbation")
    assert perturbation_feature.dtype_as_str.startswith("cat[Record[")
    assert perturbation_feature.dtype_as_object is not None

    perturbation_names = sorted(
        perturbation_feature.dtype_as_object.records.values_list("name", flat=True)
    )
    assert perturbation_names == ["DMSO", "IFNG"]

    before_count = transferred.links_feature.count()
    transferred_repeat = db.Schema.get(schema_uid).save()
    assert transferred_repeat.id == transferred.id
    assert transferred_repeat.links_feature.count() == before_count


def test_transfer_tutorial_artifact_annotations():
    db = ln.DB("laminlabs/lamindata")
    key = "example_datasets/mini_immuno/dataset1.h5ad"

    existing_local = ln.Artifact.filter(key=key).one_or_none()
    if existing_local is not None:
        existing_local.delete(storage=False, permanent=True)

    artifact = db.Artifact.get(key=key)
    artifact.save()

    artifact = db.Artifact.get(key=key)
    artifact.save(transfer="annotations")

    assert artifact.features.slots
    for schema in artifact.features.slots.values():
        # accessing schema.index should not fail after transfer
        _ = schema.index
