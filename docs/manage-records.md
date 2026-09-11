# Manage records & ontologies

```{toctree}
:maxdepth: 1
:hidden:

manage-ontologies
public-ontologies
```

```{raw} html
<iframe width="560" height="315" src="https://www.youtube.com/embed/NRzVQXJaRH8?si=Eqn4dBZyFDrbcxvm" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
```

For API examples, see {class}`~lamindb.Record`.

## Edit records like sheets

The hub offers a spreadsheet-like edit experience for {class}`~lamindb.Record` types that enforce a schema:

1. Create or select "sheet" on the left navbar under `/records`, e.g., here https://lamin.ai/laminlabs/lamindata/records

   You can organize your sheets using record types, similar to creating folders that contain different groups of related sheets. This hierarchical structure is visualized in the left navigation panel — click the `+` button to create new record types:

   <div align="center">
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/e4ueQpX1chfYETo70001.png" style="width: 40%;"/>
   </div>

2. Select the type you want.
3. Click on "Edit" to add edit rows in the sheet.
4. Each value in the sheet becomes a selectable option in dropdown menus.

   You can link columns in one sheet to another sheet or registry by specifying the appropriate data type. This creates relational connections that maintain data integrity and enable powerful cross-referencing.

   In the example shown, the `time_unit` column points to the TimeUnit sheet, which defines the permissible values for this field. This relationship ensures that only valid time units can be selected.

   <div align="center">
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/aVUC7UyIkJIxwSXv0000.png" style="width: 80%;"/>
   </div>

   Column relationships include registry links to predefined registries like TimeUnit or Treatment, ontology references that connect to biological ontologies (bionty.Organism, bionty.Tissue), cross-sheet references that link to records in other custom sheets, and enum constraints that restrict values to predefined lists.

   <div align="center">
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/uE0ax5NAXqjTfgGC0000.png" style="width: 80%;"/>
   </div>

Sheets display field definitions in the headers, showing both field names and their corresponding data types:

<div align="center">
  <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/DQLooUQbCCRxHRlx0000.png" style="width: 50%;"/>
</div>

<div align="center">
  <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/KeRnttLkyhHYXVag0000.png" style="width: 80%;"/>
</div>

## Import a csv into a sheet

Click the following button in the upper right corner of a sheet page:

<div align="center">
  <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/g6SXhxeBh7mooqsv0000.png" style="width: 80%;"/>
</div>

:::{tip}
Include a `__lamindb_record_name__` column to set the name of each record (mapped to the `record.name` field and not shown as a data column).
:::

## Export records as artifacts

Click on "Export to Artifact" to save all records under a type as an artifact.

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
