# Where to store external links and IDs?

When registering data in LaminDB, you might want to store a reference link or ID to indicate the source of the collection.

We have `reference` and `reference_type` fields for this purpose, they are available for {class}`~lamindb.Collection`, {class}`~lamindb.Transform`, {class}`~lamindb.Run` and {class}`~lamindb.Record`.

```python
# !pip install lamindb
!lamin init --storage testreference
```

```python
import lamindb as ln
```

Let's say we have a few donor samples that came form Vendor X, in order to chase back the orders, I'd like to keep track the donor ids provided by the vendor:

```python
ln.Record(
    name="donor 001", reference="VX984545", reference_type="Donor ID from Vendor X"
)
```

```python
!lamin delete --force testreference
```
