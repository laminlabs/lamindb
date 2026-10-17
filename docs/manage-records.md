---
execute_via: python
---

# Manage records & ontologies

Many use cases involving records management are primarily UI-based. However, the API allows you to perform all UI actions, too, see {class}`~lamindb.Record`.
Use cases involving ontology management are largely API-based and documented here: {doc}`/manage-ontologies`.

The following video provides an overview of the interplay of flexible records and ontology management on the UI:

```{toctree}
:maxdepth: 1
:hidden:

manage-ontologies
public-ontologies
```

```{raw} html
<div align="center">
<iframe width="560" height="315" src="https://www.youtube.com/embed/NRzVQXJaRH8?si=Eqn4dBZyFDrbcxvm" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>
```

If you haven't connected a database, create one in a development directory:

```bash
lamin init
```

## Create a data record

A data record stores values for {class}`~lamindb.Feature` objects.

```python
import lamindb as ln

gc_content = ln.Feature(name="gc_content", dtype=float).save()
sample1 = ln.Record(name="Sample 1", features={"gc_content": 0.5}).save()
sample1.describe()
```

## Organize records by type

To create and use a record type call:

```python
experiment_type = ln.Record(name="Experiments", is_type=True).save()
experiment1 = ln.Record(name="Experiment 1", type=experiment_type).save()
```

You can define a feature whose `dtype` is a record type:

```python
experiment = ln.Feature(name="experiment", dtype=experiment_type).save()
```

{meth}`~lamindb.Record.to_dataframe` exports:

```python
experiment_type.to_dataframe()
```

## Constrain records with a schema

A {class}`~lamindb.Schema` turns the type into a frame and validates the features of its data records:

```python
schema = ln.Schema(
    [experiment, gc_content.with_config(optional=True)],
).save()
samplesheet = ln.Record(name="Samples", is_type=True, schema=schema).save()
```

`experiment` is required, so setting only `gc_content` raises {class}`~lamindb.errors.ValidationError`:

```python
try:
    sample2 = ln.Record(
        name="Sample 2",
        type=samplesheet,
        features={gc_content: 0.6},
    ).save()
except ln.errors.ValidationError as error:
    print(error)
```

Also moving `sample` into `samplesheet` raises since it only has `gc_content`:

```
try:
    sample1.type = samplesheet
except ln.errors.ValidationError as error:
    print(error)

# set required featured on sample1
sample1.features.set_values({
    experiment: experiment1,
    gc_content: 0.5,
})
# move sample1 into the samplesheet
sample1.type = samplesheet
sample1.save()
```

## Use a feature as an index

You can set an index for a schema via {attr}`~lamindb.Schema.index`.

```python
sample_abbr = ln.Feature(name="sample_abbr", dtype=str).save()
score = ln.Feature(name="score", dtype=float).save()
indexed_schema = ln.Schema(features=[score], index=sample_abbr).save()
indexed_frame = ln.Record(
    name="Indexed samples", is_type=True, schema=indexed_schema
).save()

record = ln.Record(
    type=indexed_frame, features={"sample_abbr": "S-001", "score": 1.5}
).save()
record.name
```

Inspect the dataframe:

```python
indexed_frame.to_dataframe()
```

## Import records from a dataframe

{meth}`~lamindb.Record.from_dataframe` creates one data record per row. A new type name creates a record frame and infers its schema:

```python
import pandas as pd

records = ln.Record.from_dataframe(
    pd.DataFrame({"gc_content": [0.1, 0.2]}),
    type="my_df",
).save()
records
```

Pass an existing record frame as `type` to import into it. With no schema index, a `__lamindb_record_name__` column sets the record name:

```python
ln.Record.from_dataframe(
    pd.DataFrame(
        {
            "__lamindb_record_name__": ["Imported sample"],
            "experiment": ["Experiment 1"],
            "gc_content": [0.45],
        }
    ),
    type=samplesheet,
).save()
```

## Query records

Filter with a feature expression, or by a field such as `type`:

```python
ln.Record.filter(gc_content == 0.5).to_dataframe(include="features")
```

```python
ln.Record.filter(gc_content > 0.4).to_dataframe(include="features")
```

```python
ln.Record.filter(type=samplesheet).to_dataframe(include="features")
```

## Link transforms to a record

Click on the "Link" button in the "Transforms" card of the record page, e.g., here: https://lamin.ai/laminlabs/lamindata/record/mPiUOc67hQAw4Pgz

   <div align="center">
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/TPoeyphpn6evhJW60000.png" style="width: 80%;"/>
   </div>

Select a transform:

   <div align="center">
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/mfDYB2FrwsVpWy7T0000.png" style="width: 50%;"/>
   </div>

You can now launch the transform via the small "Launch" button, directly from the record:

   <div align="center">
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/tGS4pGEXqkXVscP90000.png" style="width: 50%;"/>
   </div>

## Export record frames as artifacts

Click on "Export to Artifact" to save all data records under a type as an artifact.

In the API, call {meth}`~lamindb.Record.to_artifact()` on a record frame:

```python
experiment_type.to_artifact()
```
