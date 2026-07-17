---
execute_via: python
---

# Will data & metadata stay in sync?

Yes. Transactions within Python across data & metadata are [ACID](https://en.wikipedia.org/wiki/ACID).

Because LaminDB manages a hybrid architecture—a relational database for metadata and object storage (like AWS S3) for actual file data—it uses a **two-phase commit** mechanism to guarantee that the database and storage never permanently diverge.

## How LaminDB achieves ACID

- **Atomicity:** Under normal operation, a file upload and its corresponding metadata registration succeed or fail together. If a standard error occurs (e.g., network failure or permission denied), the operation is cleanly rolled back, leaving no trace in the database or storage.
- **Consistency:** The storage state and the database state are never allowed to permanently diverge, even in the event of a hard crash (e.g., power loss). If the process is externally killed mid-upload, an internal flag ensures the database knows the file is incomplete, preventing users from querying corrupted data.
- **Isolation:** Metadata transactions inherit the strong isolation levels of the underlying relational database (e.g., PostgreSQL). Concurrent writers are safely managed by the database's native concurrency controls.
- **Durability:** Once a transaction is committed, the metadata is persisted in the relational DB and the file is safely written to storage.

### The Two-Phase Commit mechanism

To prevent dangling metadata or orphaned files, LaminDB executes a save operation in three steps:

1. **Phase 1:** LaminDB writes the metadata to the database with an internal flag `_storage_ongoing = True`.
2. **Phase 2:** The actual file bytes are uploaded to storage.
3. **Phase 3:** LaminDB updates the database to set `_storage_ongoing = False`.

If an upload process is externally killed during Phase 2, the artifact remains flagged with `_storage_ongoing = True`. This is visible in the UI, and the system knows the data is incomplete. You can then re-run `lamin save` or `artifact.save()` to attempt uploading the artifact a second time.

### ACID guarantees for collections

In tabular lakehouse formats like Iceberg or DuckLake, appending rows to a table means writing new Parquet files to storage and atomically updating a metadata file (or relational database, in DuckLake's case) to point to a new snapshot.

LaminDB uses an identical architectural pattern for dataset appends, but generalized beyond tabular data. When you call `Collection.append()`, LaminDB does not mutate existing files. Instead, it provides ACID guarantees and snapshot isolation through versioning:

1. **Artifact creation:** The new data is saved as a new `Artifact` (following the two-phase commit described above). Just like Iceberg writing a new Parquet file, the data is staged in storage.
2. **Collection versioning:** LaminDB creates a _new version_ of the `Collection` that links to the new artifact alongside the existing ones. This is the equivalent of Iceberg atomically updating its manifest to point to a new snapshot.

This architecture guarantees:

- **Atomicity:** The new collection version is only created if the new artifact is successfully stored and registered.
- **Isolation (Snapshot/Time Travel):** Concurrent readers querying the original collection version are completely unaffected by the append. The previous state of the collection remains addressable via its original version, providing the exact same "time travel" capabilities as Iceberg.

## Proving it in practice: Simulating failures

Here, we walk through different errors that can occur while saving artifacts & metadata records, and show that the LaminDB instance does not get corrupted.

```bash
lamin init --storage ./test-acid
```

```python
import pytest
import lamindb as ln
from upath import UPath

ln.settings.verbosity = "debug"

open("sample.fasta", "w").write(">seq1\nACGT\n")
```

### Simulating a failed upload within Python

Let's try to save an artifact to a storage location without permission.

```python
artifact = ln.Artifact("sample.fasta", key="sample.fasta")
```

To simulate an unauthorized storage location, we temporarily override the default storage root:

```python
ln.settings.storage._root = UPath("s3://nf-core-awsmegatests")
```

This raises an exception, and because the transaction is atomic, nothing gets saved to the database:

```python
with pytest.raises(PermissionError) as error:
    artifact.save()
print(error.exconly())
assert len(ln.Artifact.filter()) == 0
```

### Atomicity during bulk creation

Atomicity is guaranteed at the _individual artifact_ level. If a list of data objects is passed to `ln.save()` and the upload of one of these data objects fails, the successful uploads up to that point are maintained, and a `RuntimeError` is raised.

If the failure happens before any uploads begin (e.g., due to a validation error), nothing is saved:

```python
artifacts = [artifact, "this is not a record"]

with pytest.raises(Exception) as error:
    ln.save(artifacts)
print(error.exconly())
assert len(ln.Artifact.filter()) == 0  # nothing got saved
```

### Simulating an externally aborted upload

Let's restore a proper storage location:

```python
ln.settings.storage._root = UPath("./test-acid").absolute()
```

The save operation works:

```python
artifact.save()
```

Now, we simulate an upload that was killed mid-flight by manually setting the ongoing flag and deleting the file:

```python
artifact._storage_ongoing = True
artifact.save()
artifact.path.unlink()
assert artifact._aux == {"so": 1}  # storage/upload is ongoing
```

Because the system knows the upload is incomplete, we can safely re-run the save operation to recover:

```python
artifact = ln.Artifact("sample.fasta", key="sample.fasta").save()
```

```python
assert not artifact._storage_ongoing
assert artifact._aux is None
```

```bash tags=["hide-cell"]
rm -r ./test-acid
lamin delete --force test-acid
```
