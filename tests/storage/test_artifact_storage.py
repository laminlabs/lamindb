import shutil

import anndata as ad
import lamindb as ln
import pytest
from lamindb.errors import (
    IntegrityError,
)


def test_create_from_anndata_in_existing_cloud_storage():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = (
        "s3://lamindb-test/core"  # need to register as storage location
    )
    filepath = "s3://lamindb-test/core/scrnaseq_pbmc68k_tiny.h5ad"
    artifact = ln.Artifact.from_anndata(
        filepath, description="test_create_from_anndata_cloudpath"
    )
    assert artifact.n_observations == 70
    artifact.save()
    # check that the local filepath has been cleared
    assert not hasattr(artifact, "_local_filepath")
    ln.settings.storage = previous_storage


@pytest.mark.parametrize(
    "filepath_str",
    ["s3://lamindb-ci/test-data/test.parquet", "s3://lamindb-ci/test-data/test.csv"],
)
@pytest.mark.parametrize("skip_check_exists", [False, True])
@pytest.mark.parametrize("skip_size_and_hash", [False, True])
def test_create_small_file_from_remote_path(
    filepath_str, skip_check_exists, skip_size_and_hash
):
    ln.settings.creation.artifact_skip_size_hash = skip_size_and_hash
    artifact = ln.Artifact(
        filepath_str,
        skip_check_exists=skip_check_exists,
    )
    artifact.save()
    # test cache()
    file_from_local = ln.Artifact(artifact.cache(), description="test")
    # test hash equivalency when computed on local machine
    if not skip_size_and_hash:
        assert file_from_local.hash == artifact.hash
        assert file_from_local._hash_type == "md5"
        assert artifact._hash_type == "md5"
    assert artifact.path.as_posix() == filepath_str
    assert artifact.load().iloc[0].tolist() == [
        0,
        "Abingdon island giant tortoise",
        "Chelonoidis abingdonii",
        106734,
        "ASM359739v1",
        "GCA_003597395.1",
        "Full genebuild",
        "-",
        "-",
    ]
    artifact.delete(permanent=True, storage=False)
    ln.settings.creation.artifact_skip_size_hash = False


def test_create_big_file_from_remote_path():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test/core"
    filepath_str = "s3://lamindb-test/core/human_immune.h5ad"
    artifact = ln.Artifact(filepath_str)
    assert artifact.key == "human_immune.h5ad"
    assert artifact._hash_type == "md5-2"
    ln.settings.storage = previous_storage


def test_delete_artifact_from_non_managed_storage():
    artifact = ln.Artifact(
        "s3://lamindb-dev-datasets/file-to-test-for-delete.csv",
        description="My test file to delete from non-default storage",
    ).save()
    assert artifact.storage.instance_uid != ln.setup.settings.instance.uid
    assert artifact.key is not None
    filepath = artifact.path
    with pytest.raises(IntegrityError) as e:
        artifact.delete()
    assert e.exconly().startswith(
        "lamindb.errors.IntegrityError: Cannot simply delete artifacts"
    )
    artifact.delete(storage=False, permanent=True)
    assert (
        ln.Artifact.filter(
            description="My test file to delete from non-default storage",
            branch_id=None,
        ).first()
        is None
    )
    assert filepath.exists()


def test_huggingface_paths():
    artifact_adata = ln.Artifact(
        "hf://datasets/Koncopd/lamindb-test@main/anndata/pbmc68k_test.h5ad",
        description="hf adata",
    )
    artifact_adata.save()
    assert artifact_adata.hash is not None
    assert isinstance(artifact_adata.load(), ad.AnnData)
    assert artifact_adata._cache_path.exists()
    artifact_adata._cache_path.unlink()

    artifact_pq = ln.Artifact(
        "hf://datasets/Koncopd/lamindb-test/sharded_parquet", description="hf parquet"
    )
    artifact_pq.save()
    assert artifact_pq.hash is not None
    assert len(artifact_pq.open().files) == 11
    assert artifact_pq.cache().is_dir()
    shutil.rmtree(artifact_pq._cache_path)

    artifact_adata.delete(permanent=True, storage=False)
    artifact_pq.delete(permanent=True, storage=False)


def test_gcp_paths():
    artifact_folder = ln.Artifact(
        "gs://rxrx1-europe-west4/images/test/HEPG2-08", description="Test GCP folder"
    ).save()
    assert artifact_folder.hash == "6r5Hkce0UTy7X6gLeaqzBA"
    assert artifact_folder.n_files == 14772

    artifact_file = ln.Artifact(
        "gs://rxrx1-europe-west4/images/test/HEPG2-08/Plate1/B02_s1_w1.png",
        description="Test GCP file",
    ).save()
    assert artifact_file.hash == "foEgLjmuUHO62CazxN97rA"
    cache_path = artifact_file.cache()
    assert cache_path.is_file()

    cache_path.unlink()
    artifact_folder.delete(permanent=True, storage=False)
    artifact_file.delete(permanent=True, storage=False)


def test_http_paths():
    http_path = ln.UPath(
        "https://raw.githubusercontent.com/laminlabs/lamindb/refs/heads/main/README.md"
    )
    artifact_readme = ln.Artifact(http_path, description="register http readme").save()
    # might change
    assert artifact_readme.hash is not None
    cache_path = artifact_readme.cache()
    assert cache_path.exists()
    assert cache_path.stat().st_size == http_path.stat().st_size
    cache_path.unlink()
    # just check saving for the second time (when Strage record is in the db)
    artifact_license = ln.Artifact(
        "https://raw.githubusercontent.com/laminlabs/lamindb/refs/heads/main/LICENSE",
        description="register http license",
    ).save()
    assert artifact_license.hash == "IQxRSNjvb7w2OLFeWqYlsg"

    artifact_readme.delete(permanent=True, storage=False)
    artifact_license.delete(permanent=True, storage=False)


# also see test in lamindb-setup/tests/storage/test_storage_stats.py
# there is also a test for GCP there
def test_folder_like_artifact_s3():
    study0_data = ln.Artifact("s3://lamindata/iris_studies/study0_raw_images")
    assert study0_data.hash == "IVKGMfNwi8zKvnpaD_gG7w"
    assert study0_data._hash_type == "md5-d"
    assert study0_data.n_files == 51
    assert study0_data.size == 658465
