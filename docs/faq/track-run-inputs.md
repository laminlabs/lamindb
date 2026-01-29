---
executable_via: python
---

# Can I disable tracking run inputs?

Yes, if you switch {attr}`~lamindb.core.Settings.track_run_inputs` to `False`.

```python
# pip install lamindb
!lamin init --storage test-run-inputs
```

```python
import lamindb as ln
```

Some test artifacts:

```python
ln.track(transform=ln.Transform(key="Dummpy pipeline"))
ln.Artifact(ln.examples.datasets.file_jpg_paradisi05(), description="My image").save()
ln.Artifact(ln.examples.datasets.file_mini_csv(), description="My csv").save()
```

Call `ln.track()`:

```python
ln.track("Rx2s9aPTMQLY0000")
```

## Don't track artifact as run input

```python
ln.settings.track_run_inputs = False
```

```python
artifact = ln.Artifact.get(description="My image")
```

```python
artifact.cache()
```

No run inputs are linked to the current notebook run:

```python
ln.Run.get(id=ln.context.run.id).input_artifacts.all()
```

```python
artifact.view_lineage()
```

```python
assert len(ln.Run.get(id=ln.context.run.id).input_artifacts.all()) == 0
```

## Manually track artifact as run input

Let us manually track an artifact by passing `is_run_input` to either `.cache()`, `.load()` or `.open()`:

```python
artifact.cache(is_run_input=True)
```

You can see the fcs artifact is now being added to the run inputs:

```python
for input in ln.Run.get(id=ln.context.run.id).input_artifacts.all():
    print(input)
```

```python
artifact.view_lineage()
```

```python
assert len(ln.Run.get(id=ln.context.run.id).input_artifacts.all()) == 1
```

## Automatically track artifacts as run input

If you switch the following setting, and call to `.load()`, `.cache()` and `.open()` will track the artifact as run input.

```python
ln.settings.track_run_inputs = True
```

```python
artifact = ln.Artifact.get(description="My csv")
```

```python
artifact.load()
```

```python
for input in ln.Run.get(id=ln.context.run.id).input_artifacts.all():
    print(input)
```

```python
artifact.view_lineage()
```

```python
assert len(ln.Run.get(id=ln.context.run.id).input_artifacts.all()) == 2
```

```python
!lamin delete --force test-run-inputs
```
