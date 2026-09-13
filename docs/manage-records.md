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

## Export sheets as artifacts

Click on "Export to Artifact" to save all records under a type as an artifact.

In the API, call {meth}`~lamindb.Record.to_artifact()` for a `sheet`.
