# Pydantic & Pandera vs. LaminDB

This doc explains conceptual differences between data validation with `pydantic`, `pandera`, and `LaminDB`.

```python
!lamin init --storage test-pydantic-pandera --modules bionty
```

Let us work with a test dataframe.

```python
import pandas as pd
import pydantic
import lamindb as ln
import bionty as bt
import pandera.pandas as pandera
import pprint

from typing import Literal, Any

df = ln.examples.datasets.mini_immuno.get_dataset1()
df
```

## Define a schema

### pydantic

```python
Perturbation = Literal["DMSO", "IFNG"]
CellType = Literal["T cell", "B cell"]
OntologyID = Literal["EFO:0008913"]


class ImmunoSchema(pydantic.BaseModel):
    perturbation: Perturbation
    cell_type_by_model: CellType
    cell_type_by_expert: CellType
    assay_oid: OntologyID
    concentration: str
    treatment_time_h: int
    donor: str | None

    class Config:
        title = "My immuno schema"
```

### pandera

```python
pandera_schema = pandera.DataFrameSchema(
    {
        "perturbation": pandera.Column(
            str, checks=pandera.Check.isin(["DMSO", "IFNG"])
        ),
        "cell_type_by_model": pandera.Column(
            str, checks=pandera.Check.isin(["T cell", "B cell"])
        ),
        "cell_type_by_expert": pandera.Column(
            str, checks=pandera.Check.isin(["T cell", "B cell"])
        ),
        "assay_oid": pandera.Column(str, checks=pandera.Check.isin(["EFO:0008913"])),
        "concentration": pandera.Column(str),
        "treatment_time_h": pandera.Column(int),
        "donor": pandera.Column(str, nullable=True),
    },
    name="My immuno schema",
)
```

### LaminDB

Features & labels are defined on the level of the database instance.
You can either define a schema with required (and optional) columns.

```python
ln.Record(name="DMSO").save()
ln.Record(name="IFNG").save()

# leverage ontologies through types ln.Record, bt.CellType, bt.ExperimentalFactor
lamindb_schema = ln.Schema(
    name="My immuno schema",
    features=[
        ln.Feature(name="perturbation", dtype=ln.Record).save(),
        ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
        ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
        ln.Feature(name="assay_oid", dtype=bt.ExperimentalFactor.ontology_id).save(),
        ln.Feature(name="concentration", dtype=str).save(),
        ln.Feature(name="treatment_time_h", dtype=int).save(),
        ln.Feature(name="donor", dtype=str, nullable=True).save(),
    ],
).save()
```

Or merely define a constraint on the feature identifier.

```python
lamindb_schema_only_itype = ln.Schema(
    name="Allow any valid features & labels", itype=ln.Feature
)
```

## Validate a dataframe

### pydantic

```python
class DataFrameValidationError(Exception):
    pass


def validate_dataframe(df: pd.DataFrame, model: type[pydantic.BaseModel]):
    errors = []

    for i, row in enumerate(df.to_dict(orient="records")):
        try:
            model(**row)
        except pydantic.ValidationError as e:
            errors.append(f"row {i} failed validation: {e}")

    if errors:
        error_message = "\n".join(errors)
        raise DataFrameValidationError(
            f"DataFrame validation failed with the following errors:\n{error_message}"
        )
```

```python
try:
    validate_dataframe(df, ImmunoSchema)
except DataFrameValidationError as e:
    print(e)
```

To fix the validation error, we need to update the `Literal` and re-run the model definition.

```python
Perturbation = Literal["DMSO", "IFNG"]
CellType = Literal[
    "T cell", "B cell", "CD8-positive, alpha-beta T cell"  # <-- updated
]
OntologyID = Literal["EFO:0008913"]


class ImmunoSchema(pydantic.BaseModel):
    perturbation: Perturbation
    cell_type_by_model: CellType
    cell_type_by_expert: CellType
    assay_oid: OntologyID
    concentration: str
    treatment_time_h: int
    donor: str | None

    class Config:
        title = "My immuno schema"
```

```python
validate_dataframe(df, ImmunoSchema)
```

### pandera

```python
try:
    pandera_schema.validate(df)
except pandera.errors.SchemaError as e:
    print(e)
```

### LaminDB

Because the term `"CD8-positive, alpha-beta T cell"` is part of the public `CellType` ontology, validation passes the first time.

If validation had not passed, we could have resolved the issue simply by adding a new term to the `CellType` registry rather than editing the code.
This also puts downstream data scientists into a position to update ontologies.

```python
curator = ln.curators.DataFrameCurator(df, lamindb_schema)
curator.validate()
```

What was the cell type validation based on? Let's inspect the `CellType` registry.

```python
bt.CellType.to_dataframe()
```

The `CellType` regsitry is hierachical as it contains the Cell Ontology.

```python
bt.CellType.get(name="CD8-positive, alpha-beta T cell").view_parents()
```

## Overview of validation properties

Importantly, LaminDB offers not only a `DataFrameCurator`, but also a `AnnDataCurator`, `MuDataCurator`, `SpatialDataCurator`, and `TiledbsomaCurator`.

The below overview only concerns validating dataframes.

### Experience of data engineer

| property                                                                                                                       | `pydantic`                                            | `pandera`                                             | `lamindb`                                                                             |
| ------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------- | ----------------------------------------------------- | ------------------------------------------------------------------------------------- |
| define schema as code                                                                                                          | yes, in form of a `pydantic.BaseModel`                | yes, in form of a `pandera.DataFrameSchema`           | yes, in form of a `lamindb.Schema`                                                    |
| define schema as a set of constraints without the need of listing fields/columns/features; e.g. useful if validating 60k genes | no                                                    | no                                                    | yes                                                                                   |
| update labels independent of code                                                                                              | not possible because labels are enums/literals        | not possible because labels are hard-coded in `Check` | possible by adding new terms to a registry                                            |
| built-in validation from public ontologies                                                                                     | no                                                    | no                                                    | yes                                                                                   |
| sync labels with ELN/LIMS registries without code change                                                                       | no                                                    | no                                                    | yes                                                                                   |
| can re-use fields/columns/features across schemas                                                                              | limited via subclass                                  | only in same Python session                           | yes because persisted in database                                                     |
| schema modifications can invalidate previously validated datasets                                                              | yes                                                   | yes                                                   | no because LaminDB allows to query datasets that were validated with a schema version |
| can use columnar organization of dataframe                                                                                     | no, need to iterate over potentially millions of rows | yes                                                   | yes                                                                                   |

### Experience of data consumer

| property                                    | `pydantic`                                                                    | `pandera`             | `lamindb`                              |
| ------------------------------------------- | ----------------------------------------------------------------------------- | --------------------- | -------------------------------------- |
| dataset is queryable / findable             | no                                                                            | no                    | yes, by querying for labels & features |
| dataset is annotated                        | no                                                                            | no                    | yes                                    |
| user knows what validation constraints were | no, because might not have access to code and doesn't know which code was run | no (same as pydantic) | yes, via `artifact.schema`             |

## Annotation & queryability

### Engineer: annotate the dataset

Either use the `Curator` object:

```python
artifact = curator.save_artifact(key="our_datasets/dataset1.parquet")
```

If you don't expect a need for Curator functionality for updating ontologies and standardization, you can also use the `Artifact` constructor.

```python
artifact = ln.Artifact.from_dataframe(
    df, key="our_datasets/dataset1.parquet", schema=lamindb_schema
).save()
```

### Consumer: see annotations

```python
artifact.describe()
```

### Consumer: query the dataset

```python
ln.Artifact.filter(perturbation="IFNG").to_dataframe()
```

### Consumer: understand validation

By accessing `artifact.schema`, the consumer can understand _how_ the dataset was validated.

```python
artifact.schema
```

```python
artifact.schema.features.to_dataframe()
```

## Nested data with dynamic keys

We will now examine another more complex example where data is nested with potentially arbitrary (dynamic) keys.
The example is inspired by the [CELLxGENE schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/6.0.0/schema.md#uns-dataset-metadata) where annotations are stored as dictionaries in the AnnData `.uns` slot.

```python
uns_dict = ln.examples.datasets.dict_cellxgene_uns()
pprint.pprint(uns_dict)
```

### pydantic

Pydantic is primed to deal with nested data.

```python
class Images(pydantic.BaseModel):
    fullres: str
    hires: str


class Scalefactors(pydantic.BaseModel):
    spot_diameter_fullres: float
    tissue_hires_scalef: float


class Library(pydantic.BaseModel):
    images: Images
    scalefactors: Scalefactors


class Spatial(pydantic.BaseModel):
    is_single: bool
    model_config = {"extra": "allow"}

    def __init__(self, **data):
        libraries = {}
        other_fields = {}

        # store all libraries under a single key for validation
        for key, value in data.items():
            if key.startswith("library_"):
                libraries[key] = Library(**value)
            else:
                other_fields[key] = value

        other_fields["libraries"] = libraries
        super().__init__(**other_fields)


class SpatialDataSchema(pydantic.BaseModel):
    organism_ontology_term_id: str
    spatial: Spatial


validated_data = SpatialDataSchema(**uns_dict)
```

However, pydantic either requires all dictionary keys to be known beforehand to construct the Model classes or workarounds to collect all keys for a single model.

### pandera

Pandera cannot validate dictionaries because it is designed for structured dataframe data.
Therefore, we need to flatten the dictionary to transform it into a DataFrame:

```python
def _flatten_dict(d: dict[Any, Any], parent_key: str = "", sep: str = "_"):
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(_flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)
```

```python
def create_dynamic_schema(flattened_data: dict[str, Any]):
    schema_dict = {
        "organism_ontology_term_id": pandera.Column(str),
        "spatial_is_single": pandera.Column(bool),
    }

    for key in flattened_data.keys():
        if key.startswith("spatial_library_") and key.endswith("_images_fullres"):
            lib_prefix = key.replace("_images_fullres", "")
            schema_dict.update(
                {
                    f"{lib_prefix}_images_fullres": pandera.Column(str),
                    f"{lib_prefix}_images_hires": pandera.Column(str),
                    f"{lib_prefix}_scalefactors_spot_diameter_fullres": pandera.Column(
                        float
                    ),
                    f"{lib_prefix}_scalefactors_tissue_hires_scalef": pandera.Column(
                        float
                    ),
                }
            )

    return pandera.DataFrameSchema(schema_dict)


flattened = _flatten_dict(uns_dict)
df = pd.DataFrame([flattened])
spatial_schema = create_dynamic_schema(flattened)
validated_df = spatial_schema.validate(df)
```

Analogously to pydantic, pandera does not have out of the box support for dynamically named keys.
Therefore, it is necessary to dynamically construct a pydantic schema.

### LaminDB

Similarly, LaminDB currently requires constructing flattened dataframes to dynamically create features for the schema, which can then be used for validation with the DataFrameCurator.
Future improvements are expected, including support for a dictionary-specific curator.

```python
def create_dynamic_schema(flattened_data: dict[str, Any]) -> ln.Schema:
    features = []

    for key, value in flattened_data.items():
        if key == "organism_ontology_term_id":
            features.append(ln.Feature(name=key, dtype=bt.Organism.ontology_id).save())
        elif isinstance(value, bool):
            features.append(ln.Feature(name=key, dtype=bool).save())
        elif isinstance(value, (int, float)):
            features.append(ln.Feature(name=key, dtype=float).save())
        else:
            features.append(ln.Feature(name=key, dtype=str).save())

    return ln.Schema(name="Spatial data schema", features=features, coerce=True).save()


flattened = _flatten_dict(uns_dict)
flattened_df = pd.DataFrame([flattened])
spatial_schema = create_dynamic_schema(flattened)
curator = ln.curators.DataFrameCurator(flattened_df, spatial_schema)
curator.validate()
```

```{note}
Curators for scverse data structures allow for the specification of schema slots that access and validate dataframes in nested dictionary attributes like `.attrs` or `.uns`.
These schema slots use colon-separated paths like `'attrs:sample'` or `'uns:spatial:images'` to target specific dataframes for validation.
```
