---
execute_via: python
---

# Manage notebooks, scripts & workflows

Here is a 15 sec video illustrating the traceability you get by tracking data lineage.

```{raw} html
<video width="500" controls>
  <source src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/Xdiikc2c1tPtHcvF0000.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
```

**Note:** To run examples, if you don't have a `lamindb` instance, create one:

```python
!lamin init --storage ./test-track
```

## Manage notebooks and scripts

Call {meth}`~lamindb.track` to save your notebook or script as a `transform` and start tracking inputs & outputs of a run.

```{eval-rst}
.. literalinclude:: scripts/run_track_and_finish.py
   :language: python
```

<!-- #region -->

You find your notebooks and scripts in the {class}`~lamindb.Transform` registry along with pipelines & functions:

```python
transform = ln.Transform.get(key="my_analyses/my_notebook.ipynb")
transform.source_code             # source code
transform.runs.to_dataframe()     # all runs in a dataframe
transform.latest_run.report       # report of latest run
transform.latest_run.environment  # environment of latest run
```

<!-- #endregion -->

<!-- #region -->

You can use the CLI to load a transform into your current (development) directory:

```bash
lamin load --key my_analyses/my_notebook.ipynb
```

<!-- #endregion -->

<!-- #region -->

Here is how you'd load the [notebook from the video](https://lamin.ai/laminlabs/lamindata/transform/F4L3oC6QsZvQ) into your local directory:

```bash
lamin load https://lamin.ai/laminlabs/lamindata/transform/F4L3oC6QsZvQ
```

<!-- #endregion -->

### Organize local development

<!-- #region -->

If no development directory is set, script & notebook keys equal their filenames.
Otherwise, they represent the relative path in the development directory.

The exception is packaged source code, whose keys have the form `pypackages/{package_name}/path/to/file.py`.

To set the development directory to your current shell development directory, run:

```bash
lamin settings set dev-dir .
```

You can see the current status by running:

```bash
lamin info
```

<!-- #endregion -->

### Use projects

You can link the entities created during a run to a project.

```python
import lamindb as ln

my_project = ln.Project(name="My project").save()  # create & save a project
ln.track(project="My project")  # pass project
open("sample.fasta", "w").write(">seq1\nACGT\n")  # create a dataset
ln.Artifact("sample.fasta", key="sample.fasta").save()  # auto-labeled by project
```

Filter entities by project, e.g., artifacts:

```python
ln.Artifact.filter(projects=my_project).to_dataframe()
```

Access entities linked to a project:

```python
my_project.artifacts.to_dataframe()
```

The same works for `my_project.transforms` or `my_project.runs`.

### Use spaces

You can write the entities created during a run into a space that you configure on LaminHub. This is particularly useful if you want to restrict access to a space. Note that this doesn't affect bionty entities who should typically be commonly accessible.

<!-- #region -->

```python
ln.track(space="Our team space")
```

<!-- #endregion -->

<!-- #region -->

(sync-code-with-git)=

### Sync code with git

To sync scripts or workflows with their correponding files in a git repo, either export an environment variable:

```shell
export LAMINDB_SYNC_GIT_REPO = <YOUR-GIT-REPO-URL>
```

Or set the following setting:

```python
ln.settings.sync_git_repo = <YOUR-GIT-REPO-URL>
```

If you work on a single project in your lamindb instance, it makes sense to set LaminDB's `dev-dir` to the root of the local git repo clone.

```bash
dbs/
  project1/
    .git/
    script1.py
    notebook1.ipynb
  ...
```

If you work on multiple projects in your lamindb instance, you can use the `dev-dir` as the local root and nest git repositories in it.

```bash
dbs/
  database1/
    repo1/
      .git/
    repo2/
      .git/
  ...
```

<!-- #endregion -->

(manage-workflows)=

## Manage workflows

Here we'll manage workflows with `lamindb`'s {func}`~lamindb.flow` and {func}`~lamindb.step` decorators, which works out-of-the-box with the majority of Python workflow managers:

| tool      | workflow decorator | step/task decorator | notes                                          |
| --------- | ------------------ | ------------------- | ---------------------------------------------- |
| `lamindb` | `@flow`            | `@step`             | inspired by `prefect`                          |
| `prefect` | `@flow`            | `@task`             | two decorators                                 |
| `redun`   | `@task` (on main)  | `@task`             | single decorator for everything                |
| `dagster` | `@job` or `@asset` | `@op` or `@asset`   | asset-centric; `@asset` is primary             |
| `flyte`   | `@workflow`        | `@task`             | also `@dynamic` for runtime DAGs               |
| `airflow` | `@dag`             | `@task`             | TaskFlow API (modern); also supports operators |
| `zenml`   | `@pipeline`        | `@step`             | inspired by `prefect`                          |

If you're looking for more in-depth examples or for integrating with non-decorator-based workflow managers such as Nextflow or Snakemake, see {doc}`docs:pipelines`.

| tool        | workflow           | step/task         | notes            |
| ----------- | ------------------ | ----------------- | ---------------- |
| `nextflow`  | `workflow` keyword | `process` keyword | groovy-based DSL |
| `snakemake` | `rule` keyword     | `rule` keyword    | file-based DSL   |
| `metaflow`  | `FlowSpec`         | `@step`           | class-based      |
| `kedro`     | `Pipeline()`       | `node()`          | function-based   |

### A one-step workflow

Decorate a function with {func}`~lamindb.flow` to track it as a workflow:

```{eval-rst}
.. literalinclude:: scripts/my_workflow.py
   :language: python
   :caption: my_workflow.py
```

Let's run the workflow:

```python
!python scripts/my_workflow.py
```

Query the workflow via its filename:

```python
transform = ln.Transform.get(key="my_workflow.py")
transform.describe()
```

The run stored the parameter value for `key`:

```python
transform.latest_run.describe()
```

It links output artifacts:

```python
transform.latest_run.output_artifacts.to_dataframe()
```

You can query for all runs that ran with that parameter:

```python
ln.Run.filter(
    params__key="my_analysis/dataset.parquet",
).to_dataframe()
```

You can also pass complex parameters and features, see: {ref}`track-run-parameters`.

### A multi-step workflow

Here, the workflow calls an additional processing step:

```{eval-rst}
.. literalinclude:: scripts/my_workflow_with_step.py
   :language: python
   :caption: my_workflow_with_step.py
```

Let's run the workflow:

```python
!python scripts/my_workflow_with_step.py
```

The lineage of the subsetted artifact resolves the subsetting step:

```python
subsetted_artifact = ln.Artifact.get(key="my_analysis/dataset_subsetted.parquet")
subsetted_artifact.view_lineage()
```

This is the run that created the subsetted_artifact:

```python
subsetted_artifact.run
```

This is the initating run that triggered the function call:

```python
subsetted_artifact.run.initiated_by_run
```

These are the parameters of the run:

```python
subsetted_artifact.run.params
```

These are the input artifacts:

```python
subsetted_artifact.run.input_artifacts.to_dataframe()
```

These are output artifacts:

```python
subsetted_artifact.run.output_artifacts.to_dataframe()
```

### A workflow with CLI arguments

Let's use `click` to parse CLI arguments:

```{eval-rst}
.. literalinclude:: scripts/my_workflow_with_click.py
   :language: python
   :caption: my_workflow_with_click.py
```

Let's run the workflow:

```python
!python scripts/my_workflow_with_click.py --key my_analysis/dataset2.parquet
```

CLI arguments are tracked and accessible via `run.cli_args`:

```python
run = ln.Run.filter(transform__key="my_workflow_with_click.py").first()
run.describe()
```

Note that it doesn't matter whether you use `click`, `argparse`, or any other CLI argument parser.

(track-run-parameters)=

## Track parameters & features

We just saw that the function decorators `@ln.flow()` and `@ln.step()` track parameter values automatically. Here is how to pass parameters to `ln.track()`:

```{eval-rst}
.. literalinclude:: scripts/run_track_with_params.py
   :language: python
   :caption: run_track_with_params.py
```

Run the script.

```python
!python scripts/run_track_with_params.py  --input-dir ./mydataset --learning-rate 0.01 --downsample
```

Query for all runs that match certain parameters:

```python
ln.Run.filter(
    params__learning_rate=0.01,
    params__preprocess_params__downsample=True,
).to_dataframe()
```

Describe & get parameters:

```python
run = ln.Run.filter(params__learning_rate=0.01).order_by("-started_at").first()
run.describe()
run.params
```

You can also access the CLI arguments used to start the run directly:

```python
run.cli_args
```

You can also track run features in analogy to artifact features.

In contrast to params, features are validated against the `Feature` registry and allow to express relationships with entities in your registries.

Let's first define labels & features.

```python
experiment_type = ln.Record(name="Experiment", is_type=True).save()
experiment_label = ln.Record(name="Experiment1", type=experiment_type).save()
ln.Feature(name="s3_folder", dtype=str).save()
ln.Feature(name="experiment", dtype=experiment_type).save()
```

```python
!python scripts/run_track_with_features_and_params.py  --s3-folder s3://my-bucket/my-folder --experiment Experiment1
```

```python
ln.Run.filter(s3_folder="s3://my-bucket/my-folder").to_dataframe()
```

Describe & get feature values.

```python
run2 = ln.Run.filter(
    s3_folder="s3://my-bucket/my-folder", experiment="Experiment1"
).last()
run2.describe()
run2.features.get_values()
```

## Manage functions in scripts and notebooks

If you want more-fined-grained data lineage tracking in a script or notebook where you called `ln.track()`, you can also use the `step()` decorator.

### In a notebook

```python
@ln.step()
def subset_dataframe(
    input_artifact_key: str,
    output_artifact_key: str,
    subset_rows: int = 2,
    subset_cols: int = 2,
) -> None:
    artifact = ln.Artifact.get(key=input_artifact_key)
    dataset = artifact.load()
    new_data = dataset.iloc[:subset_rows, :subset_cols]
    ln.Artifact.from_dataframe(new_data, key=output_artifact_key).save()
```

Prepare a test dataset:

```python
df = ln.examples.datasets.mini_immuno.get_dataset1(otype="DataFrame")
input_artifact_key = "my_analysis/dataset.parquet"
artifact = ln.Artifact.from_dataframe(df, key=input_artifact_key).save()
```

Run the function with default params:

```python
ouput_artifact_key = input_artifact_key.replace(".parquet", "_subsetted.parquet")
subset_dataframe(input_artifact_key, ouput_artifact_key, subset_rows=1)
```

Query for the output:

```python
subsetted_artifact = ln.Artifact.get(key=ouput_artifact_key)
subsetted_artifact.view_lineage()
```

Re-run the function with a different parameter:

```python
subsetted_artifact = subset_dataframe(
    input_artifact_key, ouput_artifact_key, subset_cols=3
)
subsetted_artifact = ln.Artifact.get(key=ouput_artifact_key)
subsetted_artifact.view_lineage()
```

We created a new run:

```python
subsetted_artifact.run
```

With new parameters:

```python
subsetted_artifact.run.params
```

And a new version of the output artifact:

```python
subsetted_artifact.run.output_artifacts.to_dataframe()
```

### In a script

```{eval-rst}
.. literalinclude:: scripts/run_script_with_step.py
   :language: python
   :caption: run_script_with_step.py
```

```python
!python scripts/run_script_with_step.py --subset
```

```python
ln.view()
```

## The database

See the state of the database after we ran these different examples:

```python
ln.view()
```

## Manage notebook templates

<!-- #region -->

A notebook acts like a template upon using `lamin load` to load it. Consider you run:

```bash
lamin load https://lamin.ai/account/instance/transform/Akd7gx7Y9oVO0000
```

Upon running the returned notebook, you'll automatically create a new version and be able to browse it via the version dropdown on the UI.

Additionally, you can:

- label using `ULabel` or `Record`, e.g., `transform.records.add(template_label)`
- tag with an indicative `version` string, e.g., `transform.version = "T1"; transform.save()`
<!-- #endregion -->

<!-- #region -->

:::{dropdown} Saving a notebook as an artifact

Sometimes you might want to save a notebook as an artifact. This is how you can do it:

```bash
lamin save template1.ipynb --key templates/template1.ipynb --description "Template for analysis type 1" --registry artifact
```

:::

<!-- #endregion -->

A few checks at the end of this notebook:

```python
assert run.params == {
    "input_dir": "./mydataset",
    "learning_rate": 0.01,
    "preprocess_params": {"downsample": True, "normalization": "the_good_one"},
}, run.params
assert my_project.artifacts.exists()
assert my_project.transforms.exists()
assert my_project.runs.exists()
```
