import lamindb as ln


# we need this test both in the core and the storage/cloud tests
# because the internal logic that retrieves information about other instances
# depends on whether the current instance is managed on the hub
def test_reference_storage_location(ccaplog):
    ln.Artifact("s3://lamindata/iris_studies/study0_raw_images")
    assert (
        "referenced read-only storage location at s3://lamindata, is managed by instance with uid 4XIuR0tvaiXM"
        in ccaplog.text
    )
