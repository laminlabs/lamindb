---
execute_via: python
---

# Manage records & ontologies

This guide walks you from managing notes to managing structured relational metadata. You can build your ELN & LIMS, your own ontologies, and even a knowledge graph. If you instead want to build on public biological ontologies, read: {doc}`/manage-ontologies`.

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

## Save markdown notes

You can save a markdown note as a record:

```shell
lamin save <topic>/<my-note.md>
```

To load a record back as a markdown file, run:

```shell
lamin load https://lamin.ai/laminlabs/lamindata/record/mNDJgWFrkWQVW3ox
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

```python
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

To replace selected feature values while preserving other annotations, use
{meth}`~lamindb.models.FeatureManager.update_values`:

```python
sample1.features.update_values({gc_content: 0.6})
```

## Use a feature as an index

You can set an index for a schema via {attr}`~lamindb.Schema.index`:

```python
sample_abbr = ln.Feature(name="sample_abbr", dtype=str).save()
score = ln.Feature(name="score", dtype=float).save()
schema_with_index = ln.Schema(features=[score], index=sample_abbr).save()
frame_with_index = ln.Record(
    name="Indexed samples", is_type=True, schema=schema_with_index
).save()
record = ln.Record(
    features={"sample_abbr": "S-001", "score": 1.5},
    type=frame_with_index,
).save()
```

Inspect the dataframe:

```python
frame_with_index.to_dataframe()
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

(record-type-spaces)=

## Control spaces for record types and their records

A type's own `space` determines where that type lives. Its
{attr}`~lamindb.models.SQLRecordSettings.single_space` setting controls where its
direct children may live, including both data records and child types:

- `False`: no space constraint. Records use the normal space-selection defaults.
- `True`: require the type's own space.
- A saved {class}`~lamindb.Space` object: require that exact space, which may differ
  from the type's own space.

For a **new type constructed in Python**, `single_space` follows this priority:

1. An explicit override via `child.settings.single_space` after construction.
2. The immediate parent's current setting, copied once when constructing the
   child with `is_type=True, type=parent`.
3. `False`, if there is no parent setting to copy.

To override the initial value, assign `child.settings.single_space` to `False`,
`True`, or a `Space` before saving. Existing types keep their own settings.

This is a creation-time default. A data record uses only its immediate type's
current policy; it does not look up a grandparent's setting. Each newly created
type level can copy its parent's setting in turn, but later changes to one type
do not update other types. Ordinary data records do not receive a copied setting
themselves.

### Keep type definitions and data in different spaces

For example, a type can live in `Catalog` while requiring its child sheets and
their data records to live in `Lab`. In an instance with these two spaces, using
the `gc_content` feature defined above:

```{code-block} python
catalog = ln.Space.get(name="Catalog")
lab = ln.Space.get(name="Lab")

sheet_type = ln.Record(name="SampleSheets", is_type=True, space=catalog)
sheet_type.settings.single_space = lab
sheet_type.save()

lab_schema = ln.Schema([gc_content]).save()
sheet = ln.Record(
    name="Lab samples", is_type=True, type=sheet_type, schema=lab_schema
).save()
assert sheet.space == lab
assert sheet.settings.single_space == lab

records = ln.Record.from_dataframe(
    pd.DataFrame({"gc_content": [0.1, 0.2]}), type=sheet
).save()
assert all(record.space == lab for record in records)
```

The copied policy takes precedence over `ln.context.space` or
`ln.setup.settings.space`, even if either selects `all`. The same placement rules
apply to direct `ln.Record(..., type=sheet)` construction and dataframe imports.

### Parent and child combinations

In the table below:

- The parent type lives in `Catalog`.
- The active context or configured default space is `all`.
- A new child type and a data record under that child are created without an
  explicit `space` argument.
- Any child-setting override is made **after construction and before saving**.
- `Lab` and `Other` denote two distinct saved `Space` objects. All required writes
  are permitted, and settings do not change between construction and saving.

| Parent's `single_space` | Override on child | Child's own `space` | Child's resulting `single_space` | Data record's `space` |
| ----------------------- | ----------------- | ------------------- | -------------------------------- | --------------------- |
| `False`                 | No override       | `all`               | `False`                          | `all`                 |
| `False`                 | `False`           | `all`               | `False`                          | `all`                 |
| `False`                 | `True`            | `all`               | `True`                           | `all`                 |
| `False`                 | `Other`           | `all`               | `Other`                          | `Other`               |
| `True`                  | No override       | `Catalog`           | `True`                           | `Catalog`             |
| `True`                  | `False`           | `Catalog`           | `False`                          | `all`                 |
| `True`                  | `True`            | `Catalog`           | `True`                           | `Catalog`             |
| `True`                  | `Other`           | `Catalog`           | `Other`                          | `Other`               |
| `Lab`                   | No override       | `Lab`               | `Lab`                            | `Lab`                 |
| `Lab`                   | `False`           | `Lab`               | `False`                          | `all`                 |
| `Lab`                   | `True`            | `Lab`               | `True`                           | `Lab`                 |
| `Lab`                   | `Other`           | `Lab`               | `Other`                          | `Other`               |

For example, overriding a child's policy changes where its own children may live;
it does not exempt that child from its parent's policy:

```{code-block} python
child = ln.Record(name="Independent samples", is_type=True, type=sheet_type)
child.settings.single_space = False
child.save()  # The child itself still belongs in Lab.
```

An explicit `space` selects a destination but does not bypass a policy. PostgreSQL
rejects saving a child outside the space required by its immediate parent, whether
that child is a type or a data record. SQLite applies the Python defaults but does
not provide this database enforcement. Copying a setting grants no additional
write permissions.

### Changes after creation

Copying happens once, using the supplied parent's in-memory setting. Saving the
child does not copy it again. In particular:

- Existing or reused types are unchanged. To apply a policy to an existing sheet,
  set `sheet.settings.single_space` explicitly and save the sheet.
- Changing a parent's setting does not rewrite its children or move existing
  records. Future writes must still satisfy the immediate parent's saved policy,
  so an existing record in a now-disallowed space can fail on a later update.
- `True` continues to mean the child's own space; an exact `Space` continues to
  mean that fixed destination. If a type is allowed to move, these settings can
  therefore produce different destinations for its future children.
- Copying uses only the immediate parent. An existing intermediate type with the
  setting disabled does not inherit a grandparent's policy automatically.
- Reads, transfers, and assigning `.type` after construction do not copy the
  setting. This Python constructor behavior adds no inheritance to direct REST
  or SQL writes.

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
