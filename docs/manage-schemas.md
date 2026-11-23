# Manage schemas

Schemas define the expected structure of data in LaminDB. As requirements change, versioning your schemas helps you review updates, migrate existing data, and keep definitions consistent across instances. This guide shows you how to define, review, and apply schema changes.

## General workflow

In this example, you’ll define a sample schema and update it as requirements change.

Before you start, you’ll need:

- Compatible versions of LaminDB, lamindb-setup, and Lamin CLI.
- A Git repository with a remote you can push to for the commit and apply steps.
- A LaminDB development instance you can write to. Use its slug in place of `account/dev` in the examples.

### Define a plan

A plan declares the metadata definitions you want a LaminDB instance to have, such as schemas, features, and record types. It can also describe how to migrate existing data as those definitions change.

Create the file `plans/plan_samples.py` with this content:

```python
from lamindb_setup.evolve import Feature, LaminDType, PlanContext, Schema, create_plan


@create_plan
def samples(ctx: PlanContext, exports: None) -> None:
    sample = Feature("sample_id", dtype=LaminDType(str))
    value = Feature("measurement", dtype=LaminDType(float), nullable=False)
    Schema(name="samples", features=[sample, value], index=sample)
```

Add a unique ID so LaminDB can recognize the plan each time you apply it:

```bash
lamin evolve utils add-uid plans/plan_samples.py
```

Keep the generated UID when editing the plan. Note that this file imports directly from `lamindb_setup.evolve` and not from `lamindb`. This so you can validate the plan before connecting to an instance.

### Validate and preview

Check the plan declarations offline:

```bash
lamin evolve validate plans/plan_measurements.py
```

::::{dropdown} Example validation output
```text
╭─────────────────────────── Plan summary: samples ────────────────────────────╮
│                                                                              │
│  Plans                                                                       │
│  └── samples                                                                 │
│                                                                              │
│  Objects (3)                                                                 │
│                                                                              │
│  Features (2)                                                                │
│  ├── Feature 'sample_id'  (dtype=str)                                        │
│  └── Feature 'measurement'  (dtype=float, nullable=False)                    │
│                                                                              │
│  Schemas (1)                                                                 │
│  └── Schema 'samples'  (features=2, index='sample_id')                       │
│      └── requires                                                            │
│          ├── Feature 'sample_id'                                             │
│          └── Feature 'measurement'                                           │
│                                                                              │
╰──────────────────────────────────────────────────────────────────────────────╯
```
::::

Preview changes by comparing the plan with the live instance and any existing lock file:

```bash
lamin evolve plan plans/plan_measurements.py -i account/dev
```

::::{dropdown} Example preview output for an empty instance
```text
Plan: samples

  + Feature 'sample_id'  (known after apply)
  + Feature 'measurement'  (known after apply)
  + Schema 'samples'  (known after apply)

Plan: 3 to create, 0 to update, 0 to delete, 0 unchanged
```
::::

Review the proposed changes before applying them.

### Commit and apply

Commit and push the plan so the update can be traced to a Git revision:

```bash
git add plans/
git commit -m "Define the samples schema"
git push
lamin evolve apply plans/plan_measurements.py -i account/dev
```

Review the preview and confirm the prompt. The tool applies the definitions, tracks the operation, and writes `plans/plan_measurements.account-dev.lock`. Commit the generated lock file with your plans.

To verify that the instance matches the plan after applying it, run `lamin evolve plan` with `--fail-on-diff`. The command exits with an error if it finds any pending changes:

```bash
lamin evolve plan plans/plan_measurements.py -i account/dev --fail-on-diff
```

### Update

Suppose samples can now have a missing measurement value. Allow nulls by updating the feature in the same plan file:

```python
value = Feature("measurement", dtype=LaminDType(float), nullable=True)
```

Validate and preview the change, then commit and push the plan and apply it. The feature becomes nullable and keeps its UID. Commit the updated lock file.

### Handle refusals

#### The plan is missing its UID

If you create the plan file by hand and skip the UID step, `validate` reports that the file does not define `PLAN_TRANSFORM_STEM_UID`. Add the ID, then validate again:

```bash
lamin evolve utils add-uid plans/plan_measurements.py
lamin evolve validate plans/plan_measurements.py
```

#### Type narrowing

Suppose a sheet uses the `samples` schema and some rows lack a `measurement` value. The tool refuses to make the feature non-nullable without a way to handle those rows. Set `nullable=False` and, if zero is a valid replacement, import {func}`~lamindb.setup.evolve.fill_nulls` and add this line inside the plan:

```python
ctx.migrate(value, fill_nulls(0.0))
```
If every sheet row already has a value and no artifact datasets were validated against a schema that uses the feature, preview with `--allow-data-validated-changes`. Use the same flag on `apply` if the preview accepts the change. If validated datasets exist, check their contents separately. Once you confirm they have no null measurements, you can record that finding with `ctx.migrate(value, accept_non_null(reason=...))`.

#### The instance has drifted from the lock file

After applying the nullable update, suppose someone sets `measurement.nullable` to `False` directly in the instance. The plan and lock still say `True`. Running `lamin evolve plan` reports remote drift. Restore the live setting to `True`, then preview again. To make the feature non-nullable on purpose, make the change through the plan.

#### A rename is missing a declaration

If a change looks like a rename, the tool can refuse it. To rename the `samples` schema to `specimens` while keeping its UID, assign the new schema declaration to `schema` and add `ctx.moved(schema, "samples")` inside the plan.

## Concepts

### Plans, locks, and the instance

- **Plan:** the definitions you want to manage, including dependencies and migration declarations.
- **Lock file:** the managed object UIDs and configuration recorded after the last apply, separately for each instance.
- **Live instance:** the current database state, including changes made outside the plan.

Plans may export objects through `PlanExports` for other plans to reuse. Keep related definitions together and use separate plans for parts you manage independently.

### Types in a plan

These are declaration objects from {mod}`lamindb.setup.evolve`, distinct from the database records they create or update.

| Declaration   | Purpose                                                                                                                                      |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `FeatureType` | Group related features.                                                                                                                      |
| `Feature`     | Define a value's dtype, nullability, default, or categorical restrictions.                                                                   |
| `SchemaType`  | Group related schemas.                                                                                                                       |
| `Schema`      | Define members, optional members, an index, or component schemas.                                                                            |
| `RecordType`  | Define a type for vocabulary terms or other records.                                                                                         |
| `Sheet`       | Define a record type whose rows follow a schema.                                                                                             |
| `Record`      | Seed a named record, such as a vocabulary term. Existing records require an explicit `update_record` declaration for their values to change. |

`SchemaModel` provides an alternative way to declare schemas using Python annotations. `LaminDType` and the categorical helpers describe the feature dtypes.

### Kinds of changes

| Change               | Example and handling                                                                                                       |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| Create definitions   | Declare a new feature, schema, or record type.                                                                             |
| Edit metadata        | Update a feature description in place.                                                                                     |
| Relax a constraint   | Allow nulls, make a non-index member optional, or widen an integer dtype to float.                                         |
| Tighten a constraint | Require a previously optional member or add a categorical filter. A migration or a supported data check may be needed.     |
| Rename or move       | Use `ctx.moved()` to preserve the existing object's UID.                                                                   |
| Migrate stored data  | Use `ctx.migrate()` to backfill values, convert types, or update seeded records.                                           |
| Remove definitions   | Remove a declaration. References or retained values can require a removal policy. Seeded records require explicit removal. |

See the {mod}`Python reference <lamindb.setup.evolve>` for declarations, migration helpers, and tracking settings, and the [CLI reference](cli.md#manage-schemas) for command options, repairs, introspection, and CI workflows.

To see these concepts in code, run `lamin evolve init` with a plan name. It creates a Python plan file with commented examples of declarations, field options, and dependencies between plans.
