# CLI

Manage data with LaminDB instances.

## Manage connections

### lamin connect

Connect to an instance.

Pass a slug (`account/name`) or URL (`https://lamin.ai/account/name`).

`lamin connect` switches
{attr}`~lamindb.setup.core.SetupSettings.auto_connect` to `True` so that you
auto-connect in a Python session upon importing `lamindb`.

Alternatively, you can also connect in a Python session via {func}`~lamindb.connect`.

```text
Usage: lamin connect [OPTIONS] INSTANCE

Options:
  --help  Show this message and exit.
```

### lamin info

Show info about the environment, instance, branch, space, and user.

```text
Usage: lamin info [OPTIONS]

Options:
  --schema  View database schema.
  --help    Show this message and exit.
```

### lamin init

Init an instance.

```text
Usage: lamin init [OPTIONS]

Options:
  --storage TEXT  A local or remote folder (`'s3://...'` or `'gs://...'`).
                  Defaults to current working directory.
  --name TEXT     Instance name. If not passed, it will equal the folder name
                  passed to `storage`.
  --db TEXT       Database connection URL. Defaults to `None`, which implies
                  an SQLite file in the storage location.
  --modules TEXT  Comma-separated string of schema modules.
  --help          Show this message and exit.
```

### lamin disconnect

Disconnect from an instance.

Is the opposite of connecting to an instance.

```text
Usage: lamin disconnect [OPTIONS]

Options:
  --help  Show this message and exit.
```

## Load, save, create & delete data

### lamin load

Load a file or folder into the cache or working directory.

Pass a URL, `artifact`, or `transform`. For example:

```
lamin load https://lamin.ai/account/instance/artifact/e2G7k9EVul4JbfsEYAy5
lamin load artifact --key mydatasets/mytable.parquet
lamin load artifact --uid e2G7k9EVul4JbfsEYAy5
lamin load transform --key analysis.ipynb
lamin load transform --uid Vul4JbfsEYAy5
lamin load transform --uid Vul4JbfsEYAy5 --with-env
```

```text
Usage: lamin load [OPTIONS] ENTITY

Options:
  --uid TEXT  The uid for the entity.
  --key TEXT  The key for the entity.
  --with-env  Also return the environment for a tranform.
  --help      Show this message and exit.
```

### lamin save

Save a file or folder.

Example: Given a valid project name "my_project".

```
lamin save my_table.csv --key my_tables/my_table.csv --project my_project
```

By passing a `--project` identifier, the artifact will be labeled with the corresponding project.
If you pass a `--space` or `--branch` identifier, you save the artifact in the corresponding {class}`~lamindb.Space` or on the corresponding {class}`~lamindb.Branch`.

Note: Defaults to saving `.py`, `.ipynb`, `.R`, `.Rmd`, and `.qmd` as {class}`~lamindb.Transform` and
other file types and folders as {class}`~lamindb.Artifact`. You can enforce saving a file as
an {class}`~lamindb.Artifact` by passing `--registry artifact`.

```text
Usage: lamin save [OPTIONS] PATH

Options:
  --key TEXT          The key of the artifact or transform.
  --description TEXT  A description of the artifact or transform.
  --stem-uid TEXT     The stem uid of the artifact or transform.
  --project TEXT      A valid project name or uid.
  --space TEXT        A valid space name or uid.
  --branch TEXT       A valid branch name or uid.
  --registry TEXT     Either 'artifact' or 'transform'. If not passed, chooses
                      based on path suffix.
  --help              Show this message and exit.
```

### lamin create

Create a record for an entity.

Currently only supports creating a branch.

```
lamin create branch --name my_branch
```

```text
Usage: lamin create [OPTIONS] ENTITY

Options:
  --name TEXT  A name.
  --help       Show this message and exit.
```

### lamin delete

Delete an entity.

Currently supported: `branch`, `artifact`, and `instance`.

```
lamin delete instance --slug account/name
lamin delete branch --name my_branch
```

```text
Usage: lamin delete [OPTIONS] ENTITY

Options:
  --name TEXT
  --uid TEXT
  --slug TEXT
  --force      Do not ask for confirmation (only relevant for instance).
  --help       Show this message and exit.
```

## Describe, annotate & list data

### lamin describe

Describe an artifact.

```
lamin describe --key example_datasets/mini_immuno/dataset1.h5ad
```

```text
Usage: lamin describe [OPTIONS]

Options:
  --uid TEXT  The uid for the entity.
  --key TEXT  The key for the entity.
  --help      Show this message and exit.
```

### lamin annotate

Annotate an artifact or a transform.

You can annotate with projects and valid features & values.

```
lamin annotate --key raw/sample.fastq --project "My Project"
lamin annotate --key raw/sample.fastq --features perturbation=IFNG,DMSO cell_line=HEK297
lamin annotate --key my-notebook.ipynb --project "My Project"
```

```text
Usage: lamin annotate [OPTIONS]

Options:
  --key TEXT       The key of an artifact or transform.
  --uid TEXT       The uid of an artifact or transform.
  --project TEXT   A valid project name or uid.
  --features TEXT  Feature annotations. Supports: feature=value,
                   feature=val1,val2, or feature="val1","val2"
  --registry TEXT  Either 'artifact' or 'transform'. If not passed, chooses
                   based on key suffix.
  --help           Show this message and exit.
```

### lamin list

List records for an entity.

```
lamin list branch
lamin list space
```

```text
Usage: lamin list [OPTIONS] ENTITY

Options:
  --name TEXT  A name.
  --help       Show this message and exit.
```

## Configure

### lamin switch

Switch between branches or spaces.

Currently only supports creating a branch.

```text
Usage: lamin switch [OPTIONS]

Options:
  --branch TEXT  A valid branch name or uid.
  --space TEXT   A valid branch name or uid.
  --help         Show this message and exit.
```

### lamin cache

Manage cache.

```text
Usage: lamin cache [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  clear  Clear the cache directory.
  get    Get the cache directory.
  set    Set the cache directory.
```

### lamin settings

Manage settings.

Call without subcommands and options to show settings.

```text
Usage: lamin settings [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  set  Update settings.
```

### lamin migrate

Manage database schema migrations.

```text
Usage: lamin migrate [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  create  Create a new migration.
  deploy  Deploy migrations.
  squash  Squash migrations.
```

## Auth

### lamin login

Log into LaminHub.

`lamin login` prompts for your API key unless you set it via environment variable `LAMIN_API_KEY`.

You can create your API key in your account settings on LaminHub (top right corner).

After authenticating once, you can re-authenticate and switch between accounts via `lamin login myhandle`.

```text
Usage: lamin login [OPTIONS] [USER]

Options:
  --key TEXT  The legacy API key.
  --help      Show this message and exit.
```

### lamin logout

Log out of LaminHub.

```text
Usage: lamin logout [OPTIONS]

Options:
  --help  Show this message and exit.
```

## Other

### lamin get

Query metadata about an entity.

Currently still equivalent to `lamin describe`.

```text
Usage: lamin get [OPTIONS] [ENTITY]

Options:
  --uid TEXT  The uid for the entity.
  --key TEXT  The key for the entity.
  --help      Show this message and exit.
```

### lamin run

Run a compute job in the cloud.

This is an EXPERIMENTAL feature that enables to run a script on Modal.

Example: Given a valid project name "my_project".

```
lamin run my_script.py --project my_project
```

```text
Usage: lamin run [OPTIONS] FILEPATH

Options:
  --project TEXT    A valid project name or uid. When running on Modal,
                    creates an app with the same name.  [required]
  --image-url TEXT  A URL to the base docker image to use.
  --packages TEXT   A comma-separated list of additional packages to install.
  --cpu FLOAT       Configuration for the CPU.
  --gpu TEXT        The type of GPU to use (only compatible with cuda images).
  --help            Show this message and exit.
```
