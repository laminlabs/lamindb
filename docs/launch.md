# Launch workflows & pipelines

If you use LaminHub, you can run functions and pipelines directly in the cloud.

## Python workflows

<!-- a very similar section exists in track.md -->

In function signatures, you can use both Python and string serializations of data types to configure a Launch UI on the hub:

```python
@ln.flow()
def my_func(
    learning_rate: float,
    run_name: str,
    started_at: datetime,
    organism: "cat[bionty.Organism[source__uid=4eeXrDKBKo]]",
    sheet: "cat[Record[YSS8VU4eeXrDKBKo, is_type=True, schema__uid=6pjoBrrz4f1EzQMO]]",
    diseases: "list[cat[bionty.Disease[source__uid=4a3ejKuf]]]",
    gene_id: "cat[bionty.Gene.ensembl_gene_id[source__uid=6w75X9zM]]",
) -> str:
    ...
```

This is how it looks:

<p align="center">
  <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/1uJhxA5HshqEV8iw0000.png" alt="Lamin dtype selectors and primitive arguments" width="340" />
  <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/6GOCZBZTZwGNRClD0000.png" alt="List, artifact, and path arguments" width="340" />
</p>

See {mod}`~lamindb.base.dtypes` for more background.

## Nextflow pipelines

The hub currently only supports launching Nextflow pipelines through its Seqera integration.

The following walkthrough assumes that you already registered a pipeline in LaminDB. If you didn't, browse to the "Setup a pipeline" section. Also ensure that you have a `Compute environment` that includes a `LAMIN_API_KEY` as an env variable. It is optional but recommended to use `nf-lamin`. For this, modify your `nextflow config` as documented in [nf-lamin](nf-lamin).

### Step 0: Navigate to a pipeline

1. Navigate to the `transforms` page, e.g., [lamin.ai/laminlabs/lamindata/transforms](https://lamin.ai/laminlabs/lamindata/transforms) and select `Pipelines` to see all pipelines
2. Select your pipeline and version, e.g. [nf-core/scrnaseq](https://lamin.ai/laminlabs/lamindata/transform/PhX5TXQhj3l6wowC)
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/tSUdXKBMpm9Z7tbv0000.png" style="width: 80%;"/>
3. Click on the `Launch` button and you will land on the `launch` page. Provide your "Seqera API token" and then select `Organization`, `Workspace`, and `Compute environment`. You can also choose to "Remember this token" if you want to avoid entering the token every time.
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/3KMXObe6NemRo41j0000.png" style="width: 80%;"/>

### Step 1: Launch a pipeline

1. Optional: enter a `Run name` if you don't want to use the auto-generated name (in form of `{transform-key}-{run-number}`)
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/qzJycfLYDF4633We0000.png" style="width: 80%;"/>
2. **Parameters** (`Form` or `JSON`)

   If you prefer, you can directly edit the JSON without going through the form.

   1. Fill out `Run metadata` if configured (see setup section below for configuration)
      1. Directly fill out the form via `+ Record`
         <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/23tmjAnGkC33Ynkp0000.png" style="width: 80%;"/>
      2. Or click `Select` to select an existing row from a record frame
         <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/8WlDJtu6lWfrUXtJ0000.png" style="width: 80%;"/>
   2. Specify a file/artifact/record frame for parameters of type **file**
      1. Click on `Select` to directly select a record frame, an artifact, or a file from the s3.
         <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EAhZBrhfYDXHDlmO0000.png" style="width: 80%;"/>
      2. For parameters configured via `Input schemas`, click the `+ Create` button to fill out a new record frame.
         <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/zpPHRNTqRHyUenAK0000.png" style="width: 80%;"/>
   3. Fill out the remaining required fields marked with **\***, e.g. `outdir` for the scrnaseq pipeline
   4. Optional: provide the rest of parameters.

### Step 2: Review and submit to Seqera

1. Click `Review and Submit` at the bottom of the page
2. Review all the information you entered (can go back to edit if needed)
   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/vSPVCzF3h85QZZOV0000.png" style="width: 80%;"/>
3. Click on `Submit` to send the launch request to Seqera
4. You will see a success message if things go well (otherwise error messages) and a run id.
5. You can check the triggered run on your Seqera platform.
6. Once the pipeline run started, `nf-lamin` automatically tracks the run status. You can view them on the runs list.

### Setup a pipeline

Register a pipeline as a transform via `lamindb`:

- To register a new pipeline or a new version of an existing pipeline (please always use `Transform.from_git()`; here is the [API doc](https://docs.lamin.ai/lamindb.transform#lamindb.Transform.from_git))

  ```python
  import lamindb as ln

  # a versioned pipeline transform that points to a commit
  ln.Transform.from_git(
      url="https://github.com/nf-core/scrnaseq",
      path="main.nf",  # Path to the main script within the repository
      key="nf-core/scrnaseq",  # Optional key for the transform
      version="4.0.0",  # Optional version tag to checkout in the repository
  ).save()

  # a sliding pipeline transform that points to a branch
  ln.Transform.from_git(
      url="https://github.com/nf-core/scrnaseq",
      path="main.nf",  # Path to the main script within the repository
      key="nf-core/scrnaseq",  # Optional key for the transform
      version="dev",  # Version tag needs to equal branch
      branch="dev",  # Branch to checkout
  ).save()
  ```

Configure the settings of a pipeline:

1. Switch on the `Show configuration` toggle and you will see the `Configuration` section.
2. Click on `Settings`, you will see a pop-up window with 3 tabs: `Input schemas`, `Run metadata` and `Environment defaults`.
   1. `Input schemas`: You can configure a schema for a parameter (e.g. `input`), which allows you to create a record frame according to the schema configurations.
      <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/F4VVXZREJxTTMlxP0000.png" style="width: 80%;"/>
   2. `Run metadata`: You can configure which record frame you want to select the metadata from.
      <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/xEoxgY2NSBEucKzw0000.png" style="width: 80%;"/>
   3. `Environment defaults` : You can configure the environment defaults for the pipeline.
      <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/5shVPlptFMqymEdm0000.png" style="width: 80%;"/>

Configure pipeline sections:

1. **Configuration**

   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/qA0n4h9t9uGwbyUl0000.png" style="width: 80%;"/>
   1. `pre-run script`, `post-run script` and `nextflow config` are automatically pulled from the selected `Compute environment`
   2. Check the `nf-lamin` plugin is specified in `nextflow config` (synced from the selected compute environment)

2. **Pipeline**

   <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/Vs0VyQsuDHZA3fKp0000.png" style="width: 80%;"/>
   1. `Pipeline URL` and `Revision` are specified in transform attributes, therefore immutable here (register a new transform if you want to run a different revision of the pipeline)
   2. `Workspace directory` is pre-filled from the compute environment
   3. Optional: select `Profiles`. Once profiles are selected, you will see parameters automatically populated based on the profiles.
   4. Optional: enter a `Run name` if you don't want to use the default name
