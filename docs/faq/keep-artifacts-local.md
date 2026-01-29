# Keep artifacts local in a cloud instance

If you want to default to keeping artifacts local in a cloud instance, enable {attr}`~lamindb.setup.core.InstanceSettings.keep_artifacts_local`.

Let us first create a cloud instance that woul store artifacts exclusively on S3.

```python
!lamin login testuser1
!lamin init --storage s3://lamindb-ci/keep-artifacts-local
```

Let's import lamindb and track the current notebook run.

```python
# pip install lamindb
import lamindb as ln

ln.track("l9lFf83aPwRc")
```

## Toggling setting "keep artifacts local"

You can checkmark the "Keep artifacts local" box on the instance settings tab.

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/6Kt20kV5sQIFyV0Q0000.png" width="400px">

Or toggle it through the following instance setting.

```python
ln.setup.settings.instance.keep_artifacts_local = True
```

## Create a local storage location

Call the following for a -- potentially pre-existing -- root path and a unique host identifier.

```python
ln.Storage(root="./our_local_storage", host="abc-institute-drive1").save()
```

Now, you have two storage locations: one in the S3 bucket, and the other locally.

```python
ln.Storage.to_dataframe()
```

You can now set it as a local default storage location.
Next time you connect to the instance, this won't be necessary and the location will be automatically detected as the local default.

```python
ln.settings.local_storage = "./our_local_storage"
```

## Use a local storage location

If you save an artifact in keep-artifacts-local mode, by default, it's stored in local storage.

```python
original_filepath = ln.examples.datasets.file_fcs()
artifact = ln.Artifact(original_filepath, key="example_datasets/file1.fcs").save()
local_path = artifact.path  # local storage path
local_path
```

You'll see the `.fcs` file named by the `uid` in your `.lamindb/` directory under `./our_local_storage/`:

```python
assert artifact.path.exists()
assert artifact.path.as_posix().startswith(ln.settings.local_storage.root.as_posix())
ln.settings.local_storage.root.view_tree()
```

## Pre-existing artifacts

Assume you already have a file in your local storage location:

```python
file_in_local_storage = ln.examples.datasets.file_bam()
file_in_local_storage.rename("./our_local_storage/output.bam")
ln.UPath("our_local_storage/").view_tree()
```

When registering an artifact for it, it remains where it is.

```python
my_existing_file = ln.Artifact("./our_local_storage/output.bam").save()
ln.UPath("our_local_storage/").view_tree()
```

The storage path of the artifact matches the pre-existing file:

```python
my_existing_file.path
```

## Switching between local storage locations

You might have several local storage locations. Here is how you can switch between them.

```python
ln.Storage(root="./our_local_storage2", host="abc-institute-drive1").save()
ln.settings.local_storage = "./our_local_storage2"  # switch to the new storage location
```

Ingest a file into the new local storage location.

```python
filepath = ln.examples.datasets.file_fastq()
artifact3 = ln.Artifact(filepath, key="example_datasets/file.fastq.gz").save()
```

Inspect where all the files are.

```python
ln.Artifact.to_dataframe(include=["storage__root", "storage__region"])
```

## Upload a local artifact to the cloud

If you'd like to upload an artifact to the cloud storage location to more easily share it or view it through web applications, you pass `upload=True` to the `save()` method.

```python
artifact.save(upload=True)
```

You now see the artifact in the S3 bucket:

```python
ln.settings.storage.root.view_tree()
```

And it's no longer present in local storage:

```python
assert artifact.path.exists()
assert not local_path.exists()
assert artifact.path.as_posix().startswith(ln.settings.storage.root.as_posix())
ln.settings.local_storage.root.view_tree()
```

## Upload directly to the cloud

You can also directly upload via `upload=True`:

```python
filepath = ln.examples.datasets.file_mini_csv()
artifact2 = ln.Artifact(filepath, key="example_datasets/mini.csv").save(upload=True)
artifact2.path
```

Now we have two files on S3:

```python
ln.Artifact.to_dataframe(include="storage__root")
```

## Update storage description

You can add a description to the storage location by using the `description` field.

```python
storage_record = ln.Storage.get(root__endswith="our_local_storage")
storage_record.description = "Our shared directory for project X"
storage_record.save()
ln.Storage.to_dataframe()
```

## Delete the test instance

Delete the artifacts:

```python
artifact.delete(permanent=True)
artifact2.delete(permanent=True)
artifact3.delete(permanent=True)
my_existing_file.delete(permanent=True, storage=False)
```

Delete the instance:

```python
ln.setup.delete("keep-artifacts-local", force=True)
```
