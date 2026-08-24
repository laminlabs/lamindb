import concurrent.futures

import lamindb as ln


# we need this test both in the core and the storage/cloud tests
# because the internal logic that retrieves information about other instances
# depends on whether the current instance is managed on the hub
def test_reference_storage_location(ccaplog):
    ln.Artifact("s3://lamindata/iris_studies/study0_raw_images")
    assert ln.Storage.get(root="s3://lamindata").instance_uid == "4XIuR0tvaiXM"
    assert (
        "referenced read-only storage location at s3://lamindata, is managed by instance with uid 4XIuR0tvaiXM"
        in ccaplog.text
    )


def test_create_storage_locations_parallel():
    root: str = "nonregistered_storage"

    def create_storage() -> str:
        ln.Storage(root=root).save()  # type: ignore
        return root

    n_parallel = 3
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as executor:
        futures = [executor.submit(create_storage) for i in range(n_parallel)]
        _ = [future.result() for future in concurrent.futures.as_completed(futures)]

    storage = ln.Storage.get(root__endswith=root)
    storage.delete()
