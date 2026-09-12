# CLI

<!-- auto-generated-docs-from-here -->

Manage data with LaminDB instances.

## Configure your environment

### connect

Set the default database for this environment or directory.

This command updates your local configuration to target the specified instance:
all subsequent CLI commands and Python/R sessions will auto-connect to this instance.

You can pass a slug (`account/name`) or URL (`https://lamin.ai/account/name`).

```
# set a default instance for the current environment
lamin connect laminlabs/cellxgene
# set a default instance for the current directory
lamin connect laminlabs/cellxgene --here
# use a URL instead of a slug
lamin connect https://lamin.ai/laminlabs/cellxgene
```

Options:

```text
lamin connect [OPTIONS] INSTANCE

Options:
  --here  Connect in the current directory without changing the global default
          instance.
  --help  Show this message and exit.
```

→ Python/R alternative: create a database object via {class}`~lamindb.DB` or set the default database of your Python/R session via {func}`~lamindb.connect`

### info

Show info about the instance, development & cache directories, branch, space, and user.

Manage settings via [lamin settings](https://docs.lamin.ai/cli#settings).

Options:

```text
lamin info [OPTIONS]

Options:
  --schema  View database schema via Django plugin.
  --help    Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.setup.settings`

### init

Initialize a LaminDB instance.

Create a new development directory for your source code and `cd` into it:

```
mkdir mydata && cd mydata
```

Initialize a local SQLite database in that directory:

```
lamin init
lamin init --modules bionty
lamin init --modules bionty,pertdb
```

Initialize a SQLite database that's hosted on S3 along with all files managed by the LaminDB instance:

```
lamin init --storage s3://my-bucket
lamin init --storage gs://my-bucket
```

Initialize a PostgresSQL database with a storage location on S3:

```
lamin init --storage s3://my-bucket --db "postgresql://user:password@host:port/database"
```

Options:

```text
lamin init [OPTIONS]

Options:
  --storage TEXT  A local or remote folder (`'s3://...'` or `'gs://...'`).
                  Defaults to `./storage`.
  --name TEXT     Instance name. If no storage location is passed, uses the
                  current working directory name (like git). Otherwise uses
                  the name of the storage location.
  --db TEXT       PostgreSQL connection URI. Defaults to `None`, which implies
                  an SQLite file in the storage location.
  --modules TEXT  Comma-separated string of schema modules.
  --help          Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.setup.init`

### disconnect

Unset the default database for this environment or directory.

- Without `--here`, it clears the global default instance.
- With `--here`, it removes the nearest local marker from the current
  directory hierarchy and unsets `dev-dir` for that instance.

For example:

```
lamin disconnect
lamin disconnect --here
```

Options:

```text
lamin disconnect [OPTIONS]

Options:
  --here  Disconnect local directory context without changing the global
          default instance.
  --help  Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.setup.disconnect`

## Save, load, create & delete

### save

Save a file or folder as an `artifact`, `transform`, or `record`.

Save a **dataset** or **model** as {class}`~lamindb.Artifact`:

```
lamin save my_table.csv --key my_tables/my_table.csv
```

Pass `--store-kwargs` as a JSON object for fine-grained upload settings (normally not needed):

```
lamin save my_table.csv --key my_tables/my_table.csv --store-kwargs '{"chunksize": 1000000}'
```

Save **source code** as {class}`~lamindb.Transform`:

```
lamin save my_script.py --key my_scripts/my_script.py
```

Save a **markdown note** as {class}`~lamindb.Record`:

```
lamin save my-topic/my-note.md  # resolves `my-topic` as a record type
```

Save a **README** for the entire database:

```
lamin save README.md
```

The `save` command defaults to saving
`.py`, `.ipynb`, `.R`, `.Rmd`, and `.qmd` files as {class}`~lamindb.Transform`
and - if ommitting `--key` - `.md` files as {class}`~lamindb.Record`.
You can enforce saving a file as an {class}`~lamindb.Artifact` by passing `--registry artifact`.

You can pass a project to `--project` to label the artifact by project.
If you pass a `--space` or `--branch` identifier, you save the artifact in the corresponding {class}`~lamindb.Space` or on the corresponding {class}`~lamindb.Branch`.

Save an **agent plan** as {class}`~lamindb.Artifact`:

```
lamin save /path/to/.cursor/plans/my_task.plan.md
lamin save /path/to/.claude/plans/my_task.md
```

```{dropdown} How are agent plans handled?

Plan files are detected by suffix `.plan.md` (Cursor) or by being under `.claude/plans/`
(Claude Code). For such paths, the `key` defaults to `.plans/<filename>`, the artifact `kind`
is set to `plan`, and the description is taken from the markdown front matter (`name:` and
`overview:`). The stored artifact contains only the body (the YAML front matter is stripped).

```

**git:** When saving scripts, files will be synced with a git repo if you set:

```
export LAMINDB_SYNC_GIT_REPO=https://github.com/org/repo
```

Also see: {ref}`sync-code-with-git`

Options:

```text
lamin save [OPTIONS] PATH

Options:
  --key TEXT                      The key of the artifact or transform.
  --description TEXT              A description of the artifact or transform.
  --kind TEXT                     Artifact kind (e.g. 'plan', 'dataset',
                                  'model'). Overrides auto-inferred kind for
                                  plan files.
  --stem-uid TEXT                 The stem uid of the artifact or transform.
  --project TEXT                  A valid project name or uid.
  --space TEXT                    A valid space name or uid.
  --branch TEXT                   A valid branch name or uid.
  --registry [artifact|transform|record]
                                  Either 'artifact', 'transform', or 'record'.
                                  If not passed, chooses based on path suffix.
  --store-kwargs TEXT             Fine-grained settings for uploads as a JSON
                                  object (normally not needed), e.g.
                                  '{"chunksize": 1000000}'.
  --help                          Show this message and exit.
```

→ Python/R alternative: {class}`~lamindb.Artifact` and {class}`~lamindb.Transform`

### load

Sync a file/folder into a local cache (artifacts) or development directory (transforms).

Pass an entity or a `--key`. For example:

```
# artifacts & transforms via --key
lamin load --key mydatasets/mytable.parquet
lamin load --key analysis.ipynb
lamin load --key myanalyses/analysis.ipynb --with-env
# notes via name and topic/type hierarchy
lamin load README.md
lamin load my-topic/my-note.md
# anything via URL
lamin load https://lamin.ai/account/instance/artifact/e2G7k9EVul4JbfsE
# anything via registry and --uid
lamin load artifact --uid e2G7k9EVul4JbfsE
lamin load transform --uid Vul4JbfsEYAy5
```

Options:

```text
lamin load [OPTIONS] [ENTITY]

Options:
  --uid TEXT  The uid for the entity.
  --key TEXT  The key for the entity.
  --with-env  Also return the environment for a tranform.
  --help      Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.Artifact.load`, no equivalent for transforms

### create

Create an object.

Currently only supports creating branches and projects.

```
lamin create branch my_branch
lamin create project my_project
```

Options:

```text
lamin create [OPTIONS] {branch|project} [NAME]

Options:
  --help  Show this message and exit.
```

→ Python/R alternative: {class}`~lamindb.Branch` and {class}`~lamindb.Project`.

### delete

Delete an object.

```
# via --key or --name
lamin delete artifact --key mydatasets/mytable.parquet
lamin delete transform --key myanalyses/analysis.ipynb
lamin delete branch --name my_branch
lamin delete project --name my_project
# via --uid
lamin delete artifact --uid e2G7k9EVul4JbfsE
lamin delete transform --uid Vul4JbfsEYAy5
# via URL
lamin delete https://lamin.ai/account/db/artifact/e2G7k9EVul4JbfsE
```

To permanently delete an object, pass `--permanent`.

To delete the entire database (will ask for confirmation):

```
lamin delete account/name
```

Options:

```text
lamin delete [OPTIONS] ENTITY

Options:
  --name TEXT
  --uid TEXT
  --key TEXT   The key for the entity (artifact, transform).
  --permanent  Permanently delete the entity where applicable, e.g., for
               artifact, transform, collection.
  --force      Do not ask for confirmation (only relevant for instance).
  --help       Show this message and exit.
```

→ Python/R alternative: {meth}`~lamindb.models.SQLRecord.delete` and {func}`~lamindb.setup.delete`

## Describe, update, annotate & list

### describe

Describe an object.

Examples:

```
# via URL
lamin describe https://lamin.ai/laminlabs/lamin-site-assets/artifact/6sofuDVvTANB0f48
lamin describe https://lamin.ai/laminlabs/lamin-site-assets/transform/uDVvTANB0f48
# via --key for artifacts
lamin describe --key example_datasets/mini_immuno/dataset1.h5ad
# via registry and one of --uid / --name / --key
lamin describe artifact --uid e2G7k9EVul4JbfsE
lamin describe transform --uid Vul4JbfsEYAy5
lamin describe run --uid 6sofuDVvTANB0f48
lamin describe record --name "Experiment 1"
lamin describe project --name "My Project"
lamin describe ulabel --name "My ULabel"
lamin describe branch  # defaults to current branch
lamin describe branch --include comments
lamin describe branch --name main
```

Options:

```text
lamin describe [OPTIONS] [ENTITY]

Options:
  --uid TEXT            The uid for the entity.
  --key TEXT            The key for the entity (artifact, transform,
                        collection).
  --name TEXT           The name for the entity (record, project, ulabel,
                        branch).
  --include [comments]  Include additional content (e.g. 'comments' for readme
                        and comment blocks).
  --help                Show this message and exit.
```

→ Python/R alternative: {meth}`~lamindb.Artifact.describe`

### annotate

Annotate an artifact, transform, or collection.

You can annotate with projects, labels, records, version tags, a readme, a comment, and, for artifacts, with features. For example,

```
# via registry and --uid for any registry
lamin annotate artifact --uid e2G7k9EVul4JbfsE --project "My Project"
lamin annotate collection --uid abc123 --version "1.0"
# via registry and --name for any registry that has a name field
lamin annotate schema --name my_schema --readme README.md
# via registry and --key for any registry that as a key field
lamin annotate collection --key my_collection --version "1.0"
# via URL for any registry
lamin annotate https://lamin.ai/account/instance/artifact/e2G7k9EVul4JbfsE --project "My Project"
lamin annotate https://lamin.ai/account/instance/schema/123456ABCDEF --readme README.md
```

Annotating artifacts and transforms works via `--key` alone:

```
lamin annotate --key raw/sample.fastq --project "My Project"
lamin annotate --key raw/sample.fastq --ulabel "My ULabel" --record "Experiment 1"
lamin annotate --key raw/sample.fastq --version "1.0"
lamin annotate --key raw/sample.fastq --features perturbation=IFNG,DMSO cell_line=HEK297
lamin annotate --key raw/sample.fastq --readme README.md  # adds a readme to the artifact
lamin annotate --key raw/sample.fastq --comment "I think we should revisit this, tomorrow, WDYT?"
lamin annotate --key my-notebook.ipynb --project "My Project"
```

Branch defaults to the current branch:

```
lamin annotate branch --readme README.md  # current branch; or --name my_branch
```

Options:

```text
lamin annotate [OPTIONS] [ENTITY]

Options:
  --key TEXT       The key of an artifact, transform, or collection.
  --uid TEXT       The uid of the entity.
  --name TEXT      The name of the entity (record, project, ulabel, branch,
                   feature, schema, space).
  --project TEXT   A valid project name or uid.
  --ulabel TEXT    A valid ulabel name or uid.
  --record TEXT    A valid record name or uid.
  --version TEXT   A version tag for the artifact, transform, or collection.
  --features TEXT  Feature annotations (artifact/transform only). Supports:
                   feature=value, feature=val1,val2, or feature="val1","val2"
  --readme PATH    Path to a README file to attach as a readme block to the
                   entity.
  --comment TEXT   Comment text to attach as a comment block to the entity.
  --help           Show this message and exit.
```

→ Python/R alternative: `artifact.features.add_values()` via {meth}`~lamindb.models.FeatureManager.add_values`, `artifact.projects.add()`, `artifact.ulabels.add()`, `artifact.records.add()`, ... via {meth}`~lamindb.models.RelatedManager.add`, and `artifact.version_tag = \"1.0\"; artifact.save()` for version tags.

### update

Update mutable fields of an entity.

Examples:

```
lamin update branch --status review                  # current branch
lamin update branch --name my_branch --status draft
lamin update artifact --key my_file.parquet --description "new description"
lamin update project --name my_project --description "updated project notes"
```

Options:

```text
lamin update [OPTIONS] {artifact|transform|collection|project|branch}

Options:
  --uid TEXT                      The uid for the entity.
  --key TEXT                      The key for the entity (artifact, transform,
                                  collection).
  --name TEXT                     The name for the entity (project, branch).
  --status [standalone|draft|review|merged|closed]
                                  Set branch status (branch only).
  --description TEXT              Set description (artifact, transform,
                                  collection, project).
  --help                          Show this message and exit.
```

### get

Get a field value or describe an object.

If no field flag is passed, this behaves like `lamin describe`.
If a field flag is passed, it reads that field from the resolved entity.

Examples:

```
lamin get branch --status                # current branch status
lamin get branch --name my_branch --status
lamin get artifact --key my_file.parquet --description
```

Options:

```text
lamin get [OPTIONS] [ENTITY]

Options:
  --uid TEXT            The uid for the entity.
  --key TEXT            The key for the entity (artifact, transform,
                        collection).
  --name TEXT           The name for the entity (record, project, ulabel,
                        branch).
  --include [comments]  Include additional content (e.g. 'comments' for readme
                        and comment blocks).
  --status              Read branch status.
  --description         Read the description field.
  --help                Show this message and exit.
```

### list

List objects.

For example:

```
lamin list branch
lamin list space
```

Options:

```text
lamin list [OPTIONS] REGISTRY

Options:
  --limit INTEGER RANGE  Maximum number of rows to display.  [default: 20;
                         x>=1]
  --help                 Show this message and exit.
```

→ Python/R alternative: {meth}`~lamindb.Branch.to_dataframe()`

## Manage changes

### switch

Switch between branches.

Python/R sessions and CLI commands use the current default branch. Switch it:

```
lamin switch my_branch  # pass a name or uid of the target branch
```

To create and switch in one step, pass `-c` or `--create`:

```
lamin switch -c my_branch
```

To annotate the current branch with a `README.md`, run:

```
lamin annotate branch --readme README.md
```

To comment on the current branch, run:

```
lamin annotate branch --comment "I think we should revisit this, tomorrow, WDYT?"
```

To switch to a target space, pass `--space`:

```
lamin switch --space my_space
```

Find more info in the {class}`~lamindb.Branch` and {class}`~lamindb.Space` documents.

Options:

```text
lamin switch [OPTIONS] [TARGET]...

Options:
  --space       Switch space instead of branch.
  -c, --create  Create branch if it does not exist.
  --help        Show this message and exit.
```

→ Python/R alternative: {attr}`~lamindb.setup.core.SetupSettings.branch` and {attr}`~lamindb.setup.core.SetupSettings.space`

### merge

Merge a branch into the current branch.

Pass the `name` or `uid` of the source branch to merge into the current branch.

All `SQLRecord` objects that have `branch_id` equal to the source branch's id
are updated to the current branch's id. Example:

```
lamin switch main  # switch to the main branch
lamin merge my_branch  # after this all objects on my_branch will be on main
```

Find more info in the {class}`~lamindb.Branch` document.

Options:

```text
lamin merge [OPTIONS] BRANCH

Options:
  --help  Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.setup.merge`

## Track agents & shell scripts

### track

Track shell script runs and agent sessions.

To track a **shell script**, add `lamin track` at the beginning of the script:

```
# my_script.sh
lamin track    # initiate a tracked shell script run
lamin load --key raw/file1.txt
# do something
lamin save processed_file1.txt --key processed/file1.txt
lamin finish   # mark the tracked run as finished
```

If you run the script, input and output artifacts will be linked:

```
sh my_script.sh
```

The `lamindb` [skill](https://github.com/laminlabs/lamin-skills) ships with the `lamindb` package at `.agents/skills/`. Ask your coding agent to copy it to wherever it reads skills from — `.claude/skills/` for Claude Code, `.agents/skills/` for GitHub Copilot — so that it automatically tracks agent sessions. It will call:

```
lamin track claude   # or: lamin track copilot
# work with the agent
lamin finish
```

:::{dropdown} `lamin track copilot` says it can't find the active session?

In VS Code, make sure **"Copilot"** is selected — not **"Local"** — in the mode picker below the chat input box. `lamin track copilot` can only see sessions that go through the "Copilot"; sessions run via "Local" aren't visible to it.

```{image} https://lamin-site-assets.s3.amazonaws.com/.lamindb/f7Nw4RNYkvlw966d0000.png
:alt: Copilot mode picker
:width: 500px
```

:::

Options:

```text
lamin track [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  claude   Start tracking a Claude Code session in LaminDB.
  copilot  Start tracking a GitHub Copilot session in LaminDB.
```

→ Python/R alternative: {func}`~lamindb.track` and {func}`~lamindb.finish` for (non-shell) scripts or notebooks

### finish

Finish a tracked session.

This can be a shell script run, a Claude Code session, or a Copilot session.

Options:

```text
lamin finish [OPTIONS]

Options:
  --help  Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.finish` for (non-shell) scripts or notebooks

## Manage settings and schema & data migrations

### settings

Manage development, cache, modules, branch, space, and mount settings.

Get or set a setting by name:

- `dev-dir` → development directory {attr}`~lamindb.setup.core.SetupSettings.dev_dir`
- `cache-dir` → cache directory {attr}`~lamindb.setup.core.SetupSettings.cache_dir`
- `modules` → environment schema modules {attr}`~lamindb.setup.core.SetupSettings.modules`
- `branch` → branch {attr}`~lamindb.setup.core.SetupSettings.branch`
- `space` → space {attr}`~lamindb.setup.core.SetupSettings.space`
- `worktree` → whether dev-dir is a worktree parent

Display via [lamin info](https://docs.lamin.ai/cli#info)

Examples:

```
# dev-dir
lamin settings dev-dir get
lamin settings dev-dir set .  # set to current directory
lamin settings dev-dir set ~/my-project
lamin settings dev-dir unset
# cache-dir
lamin settings cache-dir get
lamin settings cache-dir set /path/to/cache
lamin settings cache-dir clear
# modules
lamin settings modules get
lamin settings modules set bionty,pertdb
lamin settings modules unset
# branch
lamin settings branch get
lamin settings branch set main
# space
lamin settings space get
lamin settings space set all
# worktree
lamin settings worktree get
lamin settings worktree set true
lamin settings worktree unset
# mount
lamin settings mount storage ./mnt
lamin settings mount unset ./mnt
```

Options:

```text
lamin settings [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  cache-dir  Get, set, reset, or clear the cache directory.
  dev-dir    Get or set the development directory.
  modules    Get or set environment schema modules.
  mount      Mount storage locations read-only via an installed FUSE...
  worktree   Get or set whether dev-dir is interpreted as a worktree parent.
```

→ Python/R alternative: {attr}`~lamindb.setup.core.SetupSettings.dev_dir`, {attr}`~lamindb.setup.core.SetupSettings.cache_dir`, {attr}`~lamindb.setup.core.SetupSettings.modules`, {attr}`~lamindb.setup.core.SetupSettings.branch`, and {attr}`~lamindb.setup.core.SetupSettings.space`

### migrate

Manage database schema migrations.

Options:

```text
lamin migrate [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  create  Create a new migration.
  deploy  Deploy migrations.
  squash  Squash migrations.
```

### io

Import and export databases.

Options:

```text
lamin io [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  exportdb  Export registry tables to parquet files.
  importdb  Import registry tables from parquet files.
  snapshot  Create a SQLite snapshot of the connected instance.
```

## Auth

### login

Log into LaminHub.

`lamin login` prompts for your API key unless you set it via environment variable `LAMIN_API_KEY`.

You can create your API key in your account settings on LaminHub (top right corner).

After authenticating once, you can re-authenticate and switch between accounts via `lamin login myhandle`.

Options:

```text
lamin login [OPTIONS] [USER]

Options:
  --help  Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.setup.login`

### logout

Log out of LaminHub.

Options:

```text
lamin logout [OPTIONS]

Options:
  --help  Show this message and exit.
```

## Experimental

### run

Run a compute job in the cloud.

This is an EXPERIMENTAL feature that enables to run a script on Modal.

Example: Given a valid project name "my_project",

```
lamin run my_script.py --project my_project
```

Options:

```text
lamin run [OPTIONS] FILEPATH

Options:
  --project TEXT    A valid project name or uid. When running on Modal,
                    creates an app with the same name.  [required]
  --image-url TEXT  A URL to the base docker image to use.
  --packages TEXT   A comma-separated list of additional packages to install.
  --cpu FLOAT       Configuration for the CPU.
  --gpu TEXT        The type of GPU to use (only compatible with cuda images).
  --help            Show this message and exit.
```

→ Python/R alternative: no equivalent

### hub

Query the hub API.

This is an EXPERIMENTAL feature.

Options:

```text
lamin hub [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Metadata:
  schema      Show instance schema.
  statistics  Show size and table counts.

Reads:
  list  List objects.
  get   Get one object.

Mutations:
  insert  Insert objects.
  upsert  Insert or update objects.
  update  Update objects.
  delete  Delete objects.
```
