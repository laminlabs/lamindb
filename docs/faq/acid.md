---
executable_via: python
---

# Will data & metadata stay in sync?

Here, we walk through different errors that can occur while saving artifacts & metadata records, and show that the LaminDB instance does not get corrupted by dangling metadata or artifacts.

Transactions within Python across data & metadata are [ACID](https://en.wikipedia.org/wiki/ACID).

If an upload process is externally killed and Python cannot run clean-up operations anymore, the artifact is internally still flagged with `artifact._storage_ongoing = True`. This is visible on the UI. You can then re-run `lamin save` or `artifact.save()` to attempt uploading the artifact a second time.

```python
!lamin init --storage ./test-acid
```

```python
import pytest
import lamindb as ln
from upath import UPath

ln.settings.verbosity = "debug"
```

```python
open("sample.fasta", "w").write(">seq1\nACGT\n")
```

## Save error due to failed upload within Python

Let's try to save an artifact to a storage location without permission.

```python
artifact = ln.Artifact("sample.fasta", key="sample.fasta")
```

Because the public API only allows you to set a default storage for which you have permission, we need to hack it:

```python
ln.settings.storage._root = UPath("s3://nf-core-awsmegatests")
```

This raises an exception but nothing gets saved:

```python
with pytest.raises(PermissionError) as error:
    artifact.save()
print(error.exconly())
assert len(ln.Artifact.filter()) == 0
```

## Save error during bulk creation

```python
artifacts = [artifact, "this is not a record"]
```

This raises an exception but nothing gets saved:

```python
with pytest.raises(Exception) as error:
    ln.save(artifacts)
print(error.exconly())
assert len(ln.Artifact.filter()) == 0  # nothing got saved
```

If a list of data objects is passed to `ln.save()` and the upload of one of these data objects fails, the successful uploads are maintained and a `RuntimeError` is raised, listing the successfully uploaded data objects up until that point.

## Save error due to externally aborted upload

Back to a proper storage location:

```python
ln.settings.storage._root = UPath("./test-acid").absolute()
```

The save operation works:

```python
artifact.save()
```

Let's pretend the upload was killed.

```python
artifact._storage_ongoing = True
artifact.save()
artifact.path.unlink()
assert artifact._aux == {"so": 1}  # storage/upload is ongoing
```

We can re-run it:

```python
artifact = ln.Artifact("sample.fasta", key="sample.fasta").save()
```

```python
assert not artifact._storage_ongoing
assert artifact._aux is None
```

```python
!rm -r ./test-acid
!lamin delete --force test-acid
```
