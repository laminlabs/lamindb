import lamindb as ln


def test_transfer():
    db1 = ln.DB(f"{ln.setup.settings.user.handle}/testdb1")
    artifact = db1.Artifact.get(key="README.md")
    artifact.save()
    assert artifact.key == "README.md"
    assert artifact.run is not None
    assert artifact.run.finished_at is not None
    assert artifact.run.status == "completed"
    assert artifact.run.started_at is not None
    assert ln.setup.settings.storage.root_as_str.endswith("testdb2")
    assert artifact.storage.root.endswith("testdb1")


def test_schema_transfer_ulabel_dtype():
    user_handle = ln.setup.settings.user.handle

    ln.connect("testdb1")
    perturbation_type = ln.ULabel(name="PerturbationTransferTest", is_type=True).save()
    ln.ULabel(name="DMSO", type=perturbation_type).save()
    ln.ULabel(name="IFNG", type=perturbation_type).save()
    perturbation_feature = ln.Feature(
        name="perturbation", dtype=perturbation_type
    ).save()
    schema_uid = (
        ln.Schema(
            name="transfer_schema_ulabel_perturbation",
            features=[perturbation_feature],
        )
        .save()
        .uid
    )

    ln.connect("testdb2")
    db1 = ln.DB(f"{user_handle}/testdb1")
    remote_schema = db1.Schema.get(schema_uid)
    remote_member_names = remote_schema.members.to_list("name")
    assert len(remote_member_names) > 0

    existing_local = ln.Schema.filter(uid=schema_uid).one_or_none()
    if existing_local is not None:
        existing_local.delete(permanent=True)

    existing_local_perturbation_type = ln.ULabel.filter(
        name="PerturbationTransferTest", is_type=True
    ).one_or_none()
    if existing_local_perturbation_type is not None:
        ln.ULabel.filter(type=existing_local_perturbation_type).delete(permanent=True)
        existing_local_perturbation_type.delete(permanent=True)

    assert ln.ULabel.filter(name="PerturbationTransferTest", is_type=True).count() == 0
    assert ln.ULabel.filter(type__name="PerturbationTransferTest").count() == 0

    record_only = db1.Schema.get(schema_uid).save(transfer="record")
    assert record_only.members.count() == 0

    transferred = db1.Schema.get(schema_uid).save()
    assert transferred.members.to_list("name") == remote_member_names

    perturbation = transferred.members.get(name="perturbation")
    assert perturbation.dtype_as_str.startswith("cat[ULabel[")
    assert perturbation.dtype_as_object is not None
    assert perturbation.dtype_as_object.name == "PerturbationTransferTest"

    perturbation_labels = sorted(
        perturbation.dtype_as_object.ulabels.values_list("name", flat=True)
    )
    assert perturbation_labels == ["DMSO", "IFNG"]

    before_count = transferred.links_feature.count()
    transferred_repeat = db1.Schema.get(schema_uid).save()
    assert transferred_repeat.id == transferred.id
    assert transferred_repeat.links_feature.count() == before_count


def test_schema_transfer_feature_uid_conflict_by_name():
    user_handle = ln.setup.settings.user.handle

    ln.connect("testdb1")
    source_feature = ln.Feature(name="tissue", dtype=str).save()
    source_feature_uid = source_feature.uid
    source_schema = ln.Schema(
        name="transfer_schema_feature_uid_conflict",
        features=[source_feature],
    ).save()
    schema_uid = source_schema.uid
    assert source_schema.members.get(name="tissue").uid == source_feature_uid

    ln.connect("testdb2")
    db1 = ln.DB(f"{user_handle}/testdb1")

    existing_local = ln.Schema.filter(uid=schema_uid).one_or_none()
    if existing_local is not None:
        existing_local.delete(permanent=True)

    existing_tissue_features = ln.Feature.filter(name="tissue")
    if existing_tissue_features.exists():
        existing_tissue_features.delete(permanent=True)

    local_tissue = ln.Feature(name="tissue", dtype=str).save()
    assert local_tissue.uid != source_feature_uid

    transferred = db1.Schema.get(schema_uid).save()
    transferred_tissue = transferred.members.get(name="tissue")
    assert transferred_tissue.uid == source_feature_uid
    assert transferred_tissue.uid != local_tissue.uid
