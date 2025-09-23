# .latest_version is tested in test_versioning.py


import re

import bionty as bt
import lamindb as ln
import pytest
from django.core.exceptions import FieldError
from lamindb.base.users import current_user_id
from lamindb.errors import InvalidArgument
from lamindb.models import ArtifactSet, BasicQuerySet, QuerySet
from lamindb.models.query_set import DoesNotExist


# please also see the test_curate_df.py tests
def test_df():
    project_label = ln.ULabel(name="project").save()
    project_names = [f"Project {i}" for i in range(3)]
    labels = ln.ULabel.from_values(project_names, create=True).save()
    project_label.children.add(*labels)
    df = ln.ULabel.to_dataframe(include="parents__name")
    assert df.columns[2] == "parents__name"
    assert df["parents__name"].iloc[0] == {project_label.name}
    df = ln.ULabel.to_dataframe(include=["parents__name", "parents__created_by_id"])
    assert df.columns[3] == "parents__created_by_id"
    assert df["parents__name"].iloc[0] == {project_label.name}
    assert set(df["parents__created_by_id"].iloc[0]) == {current_user_id()}

    # for other models
    feature_names = [f"Feature {i}" for i in range(3)]
    features = [ln.Feature(name=name, dtype=int) for name in feature_names]
    ln.save(features)
    schema = ln.Schema(features, name="my schema").save()
    schema.features.set(features)

    df = ln.Schema.filter(name="my schema").to_dataframe(include="features__name")
    assert df.columns[2] == "features__name"
    # order is not conserved
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    # pass a list
    df = ln.Schema.filter(name="my schema").to_dataframe(
        include=["features__name", "features__created_by_id"]
    )
    assert df.columns[3] == "features__created_by_id"
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    assert set(df["features__created_by_id"].iloc[0]) == {current_user_id()}

    # inner join parents on features
    df = ln.Schema.filter().to_dataframe(
        include=["features__name", "features__created_by_id"]
    )
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    assert set(df["features__created_by_id"].iloc[0]) == {current_user_id()}

    # raise error for non many-to-many
    df = ln.ULabel.filter(name="Project 0").to_dataframe(include="created_by__name")
    assert df["created_by__name"].iloc[0] == "Test User1"

    # do not return fields with no data in the registry
    # does not make sense in Alex's opinion
    # too much magic; got removed in https://github.com/laminlabs/lamindb/pull/2238
    # df = (
    #     ln.Artifact.using("laminlabs/cellxgene")
    #     .filter(suffix=".h5ad")
    #     .to_dataframe(include=["tissues__name", "pathways__name"])
    # )
    # assert "tissues__name" in df.columns
    # assert "pathways__name" not in df.columns
    # assert df.shape[0] > 0

    # clean up
    project_label.delete(permanent=True)
    for label in labels:
        label.delete(permanent=True)

    schema.delete(permanent=True)
    for feature in features:
        feature.delete(permanent=True)

    # call it from a non-select-derived queryset
    qs = ln.User.objects.all()
    assert qs.to_dataframe().iloc[0]["handle"] == "testuser1"


def test_one_first():
    qs = ln.User.objects.all()
    assert qs.one().handle == "testuser1"
    assert qs.first().handle == "testuser1"
    assert qs.one_or_none().handle == "testuser1"

    qs = ln.User.filter(handle="test")
    with pytest.raises(DoesNotExist):
        qs.one()
    qs = bt.Source.filter().all()
    with pytest.raises(ln.errors.MultipleResultsFound):
        qs.one()
    with pytest.raises(ln.errors.MultipleResultsFound):
        qs.one_or_none()


def test_filter_related_field_name():
    with pytest.raises(
        FieldError,
        match=re.escape(
            "Invalid lookup 'somelabel' for ulabels. Did you mean ulabels__name?"
        ),
    ):
        ln.Artifact.filter(ulabels="somelabel").all()


def test_filter_unknown_field():
    with pytest.raises(InvalidArgument) as error:
        ln.Artifact.filter(nonexistent="value").all()
    assert error.exconly() == (
        "lamindb.errors.InvalidArgument: You can query either by available fields: branch, cell_lines, cell_markers, cell_types, collections, created_at, created_by, description, developmental_stages, diseases, ethnicities, experimental_factors, feature_sets, genes, hash, id, input_of_runs, is_latest, key, kind, linked_in_records, n_files, n_observations, organisms, otype, pathways, phenotypes, projects, proteins, records, references, run, schema, size, space, storage, suffix, tissues, transform, uid, ulabels, updated_at, version, visibility\n"
        "Or fix invalid feature names: nonexistent"
    )


def test_get_id_type_error():
    with pytest.raises(
        ValueError, match=re.escape("Field 'id' expected a number but got 'abc'.")
    ):
        ln.Artifact.get(id="abc")


def test_get_related_field_name():
    with pytest.raises(
        FieldError,
        match=re.escape(
            "Invalid lookup 'somelabel' for ulabels. Did you mean ulabels__name?"
        ),
    ):
        ln.Artifact.get(ulabels="somelabel").all()


def test_get_unknown_field():
    with pytest.raises(
        FieldError,
        match=re.escape(
            "Unknown field 'nonexistent'. Available fields: branch, cell_lines, cell_markers, cell_types, collections, created_at, created_by, description, developmental_stages, diseases, ethnicities, experimental_factors, feature_sets, genes, hash, id, input_of_runs, is_latest, key, kind, linked_in_records, n_files, n_observations, organisms, otype, pathways, phenotypes, projects, proteins, records, references, run, schema, size, space, storage, suffix, tissues, transform, uid, ulabels, updated_at, version, visibility"
        ),
    ):
        ln.Artifact.get(nonexistent="value")


def test_search():
    label_names = [f"ULabel {i}" for i in range(3)]
    labels = [ln.ULabel(name=name) for name in label_names]
    ln.save(labels)
    qs = ln.ULabel.filter(name__startswith="ULabel").all()
    assert qs.search("ULabel 1")[0].name == "ULabel 1"
    assert qs.search("ULabel 1", field=ln.ULabel.name)[0].name == "ULabel 1"
    for label in labels:
        label.delete(permanent=True)


def test_lookup():
    qs = ln.User.filter(handle="testuser1").all()
    # pass str to field
    lookup = qs.lookup(field="handle")
    assert lookup.testuser1.handle == "testuser1"
    # pass StrField to field
    lookup = qs.lookup(field=ln.User.handle)
    assert lookup.testuser1.handle == "testuser1"
    # manager, default field
    qsm = ln.User.filter(handle="testuser1")
    lookup = qsm.lookup()
    assert lookup.testuser1.handle == "testuser1"


def test_inspect():
    qs = ln.User.filter(handle="testuser1").all()
    assert qs.inspect(["user1", "user2"], "name")["validated"] == []
    assert ln.User.inspect(["user1", "user2"], "name")["validated"] == []
    assert ln.User.inspect(["user1", "user2"], ln.User.name)["validated"] == []
    assert ln.User.inspect("user1", "name")["validated"] == []


def test_validate():
    qs = ln.User.filter(handle="testuser1").all()
    assert qs.validate(["testuser1", "Test User1"], "handle").tolist() == [True, False]
    assert ln.User.validate(["testuser1", "Test User1"], "handle").tolist() == [
        True,
        False,
    ]
    assert ln.User.validate(["testuser1", "Test User1"], ln.User.handle).tolist() == [
        True,
        False,
    ]
    # returns True
    assert ln.User.validate("testuser1", ln.User.handle)


def test_map_synonyms():
    qs = ln.User.filter(handle="testuser1").all()
    assert qs.standardize(["user1", "user2"]) == ["user1", "user2"]


def test_get_doesnotexist_error():
    non_existent_label = "some-label-name"

    with pytest.raises(DoesNotExist) as excinfo:
        ln.ULabel.get(non_existent_label)

    error_message = str(excinfo.value)
    assert f"No record found with uid '{non_existent_label}'" in error_message
    assert (
        f"Did you forget a keyword as in ULabel.get(name='{non_existent_label}')?"
        in error_message
    )


def test_get():
    branch = ln.Branch(name="test_branch").save()

    artifact = ln.Artifact.from_dataframe(
        ln.User.to_dataframe(), key="df_test_get.parquet"
    )
    artifact.branch = branch
    artifact.save()
    # errors if doesn't find or multiple records found
    ln.Artifact.get(key="df_test_get.parquet", branch=branch)
    ln.Artifact.get(key="df_test_get.parquet", branch_id=branch.id)
    ln.Artifact.get(key="df_test_get.parquet", branch__in=[branch])
    ln.Artifact.get(key="df_test_get.parquet", branch_id__in=[branch.id])
    ln.Artifact.get(key="df_test_get.parquet", branch=None)
    ln.Artifact.get(key="df_test_get.parquet", branch_id=None)

    with pytest.raises(ln.Artifact.DoesNotExist):
        ln.Artifact.get(key="df_test_get.parquet")

    ln.Artifact.get(artifact.id)
    ln.Artifact.get(id=artifact.id)
    ln.Artifact.get(id__in=[artifact.id])

    ln.Artifact.get(artifact.uid[:5])
    ln.Artifact.get(uid=artifact.uid)
    ln.Artifact.get(uid__in=[artifact.uid])

    ln.Artifact.get(hash=artifact.hash)
    ln.Artifact.get(hash__in=[artifact.hash])

    artifact.delete(permanent=True)
    branch.delete()


def test_to_class():
    qs = ln.Artifact.filter()
    assert isinstance(qs, QuerySet)
    assert isinstance(qs, ArtifactSet)

    qs_copy = qs._to_non_basic(copy=True)
    assert isinstance(qs_copy, QuerySet)
    assert isinstance(qs_copy, ArtifactSet)

    qs_basic = qs._to_basic(copy=True)
    assert isinstance(qs_basic, BasicQuerySet)
    assert isinstance(qs_basic, ArtifactSet)
    assert not isinstance(qs_basic, QuerySet)

    qs_basic._to_non_basic(copy=False)
    assert isinstance(qs_basic, QuerySet)
    assert isinstance(qs_basic, ArtifactSet)
