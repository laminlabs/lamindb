# CLI

<!-- auto-generated-docs-from-here -->

Manage data with LaminDB instances.

:::{dropdown} How do I enable shell autocompletion?

On MacOS, enable shell autocompletion by adding the following to your `~/.zshrc`:

```bash
autoload -Uz compinit && compinit
eval "$(_LAMIN_COMPLETE=zsh_source lamin)"
```

:::

## Configure your environment

### connect

Set the default database for this directory.

This command updates your local configuration to target the specified instance:
all subsequent CLI commands and Python/R sessions will auto-connect to this instance.

You can pass a slug (`account/name`) or URL (`https://lamin.ai/account/name`).
Pass `--use-proxy-db` to connect through the database proxy specified in the hub.

```
# set a default database for the current directory (recommended)
lamin connect laminlabs/cellxgene --here
# set a default database for the entire home directory
lamin connect laminlabs/cellxgene
# use a URL instead of a slug
lamin connect https://lamin.ai/laminlabs/cellxgene
# connect through the database proxy specified in the hub
lamin connect laminlabs/cellxgene --use-proxy-db
```

Options:

```text
lamin connect [OPTIONS] INSTANCE

Options:
  --here          Connect in the current directory without changing the global
                  default instance.
  --use-proxy-db  Use the proxy database connection specified in the hub.
  --help          Show this message and exit.
```

→ Python/R alternative: create a database object via {class}`~lamindb.DB` or set the default database of your Python/R session via {func}`~lamindb.connect`

### info

Show info about the database, branch, space, and user.

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

Initialize a database in the current directory.

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

Unset the default database for this directory.

- With `--here`, it clears the default database for the current directory.
- Without `--here`, it clears the default database for the entire home directory.

For example:

```
lamin disconnect --here
lamin disconnect
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

Save a file or folder as an artifact, transform, record, or project.

Save a **dataset** or **model** as {class}`~lamindb.Artifact`:

```
lamin save my_table.csv --key my_tables/my_table.csv
```

Pass `--batch-size` to control parallel file transfers when saving artifact folders (default 128). Reducing to 20 or lower can help with network errors:

```
lamin save my_folder --key my_tables --batch-size 20
```

Pass `--store-kwargs` as a JSON object for fine-grained upload settings (normally not needed):

```
lamin save my_table.csv --key my_tables/my_table.csv --store-kwargs '{"chunksize": 1000000}'
```

Save **source code** as {class}`~lamindb.Transform`:

```
lamin save my_script.py --key my_scripts/my_script.py
```

Save a **markdown note** in `Projects/` as {class}`~lamindb.Project`:

```
lamin save Projects/my-project.md
lamin save Projects/my-type/my-project.md  # `my-type` is an existing project type
```

Save a **markdown note** as {class}`~lamindb.Record`:

```
lamin save my-topic/my-note.md  # resolves `my-topic` as a record page
```

Save a **README** for the entire database:

```
lamin save README.md
```

The `save` command defaults to saving
`.py`, `.ipynb`, `.R`, `.Rmd`, and `.qmd` files as {class}`~lamindb.Transform`
and - if ommitting `--key` - `.md` files as {class}`~lamindb.Record` or {class}`~lamindb.Project`.
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
  --registry [artifact|transform|record|project]
                                  Either 'artifact', 'transform', 'record', or
                                  'project'. If not passed, chooses based on
                                  path.
  --batch-size INTEGER            Number of files transferred in parallel when
                                  saving artifact folders (default 128).
                                  Reducing to 20 or lower can help with
                                  network errors.
  --store-kwargs TEXT             Fine-grained settings for uploads as a JSON
                                  object (normally not needed), e.g.
                                  '{"chunksize": 1000000}'.
  --help                          Show this message and exit.
```

→ Python/R alternative: {class}`~lamindb.Artifact` and {class}`~lamindb.Transform`

### load

Sync a file/folder into a local cache (artifacts) or development directory (transforms, records).

Pass an entity or a `--key`. For example:

```
# artifacts & transforms via --key
lamin load --key mydatasets/mytable.parquet
lamin load --key analysis.ipynb
lamin load --key myanalyses/analysis.ipynb --with-env
# notes via name and topic/type hierarchy
lamin load README.md
lamin load My-topic/my-note.md
# anything via URL
lamin load https://lamin.ai/account/instance/artifact/e2G7k9EVul4JbfsE
# anything via registry and --uid
lamin load artifact --uid e2G7k9EVul4JbfsE
lamin load transform --uid Vul4JbfsEYAy5
```

Pass `--batch-size` to control parallel file transfers when loading artifact or collection folders (default 128). Reducing to 20 or lower can help with network errors:

```
lamin load --key mydatasets/myfolder --batch-size 20
```

Pass `--store-kwargs` as a JSON object for fine-grained artifact or collection download settings (normally not needed):

```
lamin load --key mydatasets/mytable.parquet --store-kwargs '{"batch_size": 20}'
```

Options:

```text
lamin load [OPTIONS] [ENTITY]

Options:
  --uid TEXT            The uid for the entity.
  --key TEXT            The key for the entity.
  --with-env            Also return the environment for a tranform.
  --batch-size INTEGER  Number of files transferred in parallel when loading
                        artifact or collection folders (default 128). Reducing
                        to 20 or lower can help with network errors.
  --store-kwargs TEXT   Fine-grained settings for artifact or collection
                        downloads as a JSON object (normally not needed), e.g.
                        '{"batch_size": 20}'.
  --help                Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.Artifact.load`, no equivalent for transforms

### create

Create an object.

Currently only supports creating branches and projects.

Use `lamin save` to create artifacts, transforms, and records.

Examples:

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

Annotate an object.

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

Update an object.

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

Get object metadata.

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

## Manage schemas

(cli-evolve)=

### evolve

Manage schema plans, migrations, and repairs.

Options:

```text
lamin evolve [OPTIONS] COMMAND [ARGS]...

Options:
  --version            Print the installed lamindb-setup evolution version and
                       exit.
  --log-queries        Log every SQL statement issued against an instance to
                       stderr, with its duration, and a total when the command
                       finishes. Use it to see which query a slow command is
                       waiting on.
  --project-root PATH  Directory plan files are imported relative to, so a
                       plan can import another one by its full dotted path
                       ('from plans.wetlab.plan_x import y'). Defaults to the
                       nearest ancestor of each plan holding a
                       lamin_schemas.yml, pyproject.toml or .git, else the
                       plan's own directory. Pass this only when the plans are
                       scattered or vendored outside such a root.
  --profile PATH       Profile the command with pyinstrument and write the
                       report here. A '.html' suffix produces the interactive
                       flame view, anything else a plain-text tree. Needs the
                       'profile' extra.
  --help               Show this message and exit.

Commands:
  init               Scaffold a new plan file with a fresh transform stem...
  validate           Finalize and render one or more plans offline (no...
  plan               Show a read-only diff of what applying the plan(s)...
  apply              Create the plan(s)' objects in the instance and...
  destroy            Tear down the plan(s)' objects, refusing if any are...
  introspect         Generate a plan file from an existing instance.
  repair             Repair an instance: dangling objects, duplicates,...
  utils              Miscellaneous helpers for authoring plans.
  ci                 Run plan checks in CI from a lamin_schemas.yml...
  prune-dangling     Deprecated alias for 'repair dangling'. (DEPRECATED)
  remove-duplicates  Deprecated alias for 'repair duplicates'. (DEPRECATED)
  verify-dtypes      Deprecated alias for 'repair dtypes'. (DEPRECATED)
```

::::{dropdown} lamin evolve init
:name: cli-evolve-init

Scaffold a new plan file with a fresh transform stem uid.

```text
lamin evolve init [OPTIONS] [NAME]

Options:
  -d, --output-dir, --directory PATH
                                  Directory in which to create the new plan
                                  file.  [default: .]
  --add-uid PATH                  Deprecated: use 'lamin evolve utils add-uid
                                  <plan>'. Add a stem uid to this existing
                                  plan file instead of scaffolding.
  --minimal                       Scaffold a bare plan without the commented
                                  building-block examples.
  -i, --instance TEXT             Instance(s) to check the generated uid
                                  against for collisions. Repeat to require
                                  the uid to be free on several instances.
  --help                          Show this message and exit.

Arguments:
  NAME  Plan name, e.g. 'mynewplan' -> plan_mynewplan.py.
```

::::

::::{dropdown} lamin evolve validate
:name: cli-evolve-validate

Finalize and render one or more plans offline (no instance connection).

```text
lamin evolve validate [OPTIONS] TARGETS...

Options:
  -a, --all-plans                 Include every plan defined in each provided
                                  file, instead of a single plan per target.
                                  Targets must be plain file paths without a
                                  ':plan_variable' suffix.
  -f, --format [terminal|markdown|json]
                                  Output format: 'terminal' for colored output
                                  (default), 'markdown' for a GitHub-flavored
                                  fragment to post on a pull request, or
                                  'json' for a machine-readable document to
                                  pipe into a script.  [default: terminal]
  -o, --output PATH               Also write the rendered document to this
                                  file. Needs --format markdown or --format
                                  json.
  --help                          Show this message and exit.

Arguments:
  TARGETS  One or more plans to act on, each as
           'path/to/plan.py:plan_variable' (the ':plan_variable' suffix is
           optional when the file defines exactly one plan). When several are
           given they are evaluated together: an object dropped from one plan
           but still declared by another listed plan is treated as moved
           between plans, not deleted.
```

::::

::::{dropdown} lamin evolve plan
:name: cli-evolve-plan

Show a read-only diff of what applying the plan(s) would create.

```text
lamin evolve plan [OPTIONS] TARGETS...

Options:
  -i, --instance TEXT             LaminDB instance to act against, as
                                  'account/instance'.  [required]
  -f, --format [terminal|markdown|json]
                                  Output format: 'terminal' for colored output
                                  (default), 'markdown' for a GitHub-flavored
                                  fragment to post on a pull request, or
                                  'json' for a machine-readable document to
                                  pipe into a script.  [default: terminal]
  -o, --output PATH               Also write the rendered document to this
                                  file. Needs --format markdown or --format
                                  json.
  --preview-migrations            Rehearse the apply inside a transaction that
                                  is always rolled back: create what the plan
                                  creates, run every declared ctx.migrate()
                                  migration and its checks, and report the
                                  outcome in the diff -- or fail as the apply
                                  would. Nothing is persisted. Off by default
                                  so the diff stays a pure read-only preview
                                  that never runs migration code.
  --strict-uids / --no-strict-uids
                                  Judge an object whose uid differs from the
                                  uid the plan's lock files for other
                                  instances record as 'apply' does (default):
                                  --fail-on-diff counts it as pending, and
                                  --preview-migrations refuses to create an
                                  object such a lock pins but something else
                                  holds here. --no-strict-uids only reports
                                  the misalignment and rehearses the apply
                                  that creates the object with a fresh uid
                                  instead.  [default: strict-uids]
  --strict-renames / --no-strict-renames
                                  Refuse to apply a delete+create that looks
                                  like a rename not declared with ctx.moved
                                  (default). Use --no-strict-renames to allow
                                  it when the objects are genuinely distinct.
                                  'plan' shows the diff and exits 1 wherever
                                  'apply' would refuse; --no-strict-renames
                                  only warns.  [default: strict-renames]
  --allow-data-validated-changes  Allow a blocked validation-affecting change
                                  to a feature when the live data proves it
                                  safe: the feature has no stored values, or
                                  every value validates against the new dtype.
                                  Off by default. 'plan' then shows the diff
                                  'apply --allow-data-validated-changes' would
                                  apply.
  --fail-on-diff                  Exit 1 once every plan is shown (and
                                  --output written) when any of them has a
                                  pending change of its own, or has a
                                  misaligned uid 'apply' would refuse (see
                                  --strict-uids). Off by default: only a
                                  refusal fails the run.
  --show-unchanged                List every object of the plan in the
                                  terminal diff, unchanged ones included.
  --help                          Show this message and exit.

Arguments:
  TARGETS  One or more plans to act on, each as
           'path/to/plan.py:plan_variable' (the ':plan_variable' suffix is
           optional when the file defines exactly one plan). When several are
           given they are evaluated together: an object dropped from one plan
           but still declared by another listed plan is treated as moved
           between plans, not deleted.
```

::::

::::{dropdown} lamin evolve apply
:name: cli-evolve-apply

Create the plan(s)' objects in the instance and write the lock file(s).

```text
lamin evolve apply [OPTIONS] TARGETS...

Options:
  -i, --instance TEXT             LaminDB instance to act against, as
                                  'account/instance'.  [required]
  -y, --yes                       Skip the confirmation prompt.
  --strict-renames / --no-strict-renames
                                  Refuse to apply a delete+create that looks
                                  like a rename not declared with ctx.moved
                                  (default). Use --no-strict-renames to allow
                                  it when the objects are genuinely distinct.
                                  With --config, 'strict_renames: false' there
                                  allows it too.  [default: strict-renames]
  --allow-data-validated-changes  Allow a blocked validation-affecting change
                                  to a feature when the live data proves it
                                  safe: the feature has no stored values, or
                                  every value validates against the new dtype.
                                  Off by default. With --config,
                                  'allow_data_validated_changes: true' there
                                  opts in too.
  --strict-uids / --no-strict-uids
                                  Refuse to apply while an object's uid
                                  differs from the uid the plan's lock files
                                  for other instances record (default); align
                                  them with 'repair align-uids'. Use --no-
                                  strict-uids (or, with --config,
                                  'strict_uids: false' there) to apply anyway,
                                  creating an object that cannot be pinned
                                  with a fresh uid.  [default: strict-uids]
  -f, --format [terminal|markdown|json]
                                  Output format: 'terminal' for colored output
                                  (default), 'markdown' for a GitHub-flavored
                                  fragment to post on a pull request, or
                                  'json' for a machine-readable document to
                                  pipe into a script.  [default: terminal]
  -o, --output PATH               Also write the rendered document to this
                                  file. Needs --format markdown or --format
                                  json.
  -c, --config PATH               Path to a lamin_schemas.yml whose settings
                                  this apply honours as 'ci apply' does: the
                                  instance's 'track_space'/'track_prefix'
                                  (falling back to the LAMIN_SCHEMA_EVOLUTION_
                                  TRACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PRE
                                  FIX env vars),
                                  'allow_data_validated_changes' (ORed with
                                  the flag), and
                                  'strict_uids'/'strict_renames' (ANDed with
                                  the flags).
  --show-unchanged                List every object of the plan in the
                                  terminal diff, unchanged ones included.
  --help                          Show this message and exit.

Arguments:
  TARGETS  One or more plans to act on, each as
           'path/to/plan.py:plan_variable' (the ':plan_variable' suffix is
           optional when the file defines exactly one plan). When several are
           given they are evaluated together: an object dropped from one plan
           but still declared by another listed plan is treated as moved
           between plans, not deleted.
```

::::

::::{dropdown} lamin evolve destroy
:name: cli-evolve-destroy

Tear down the plan(s)' objects, refusing if any are referenced.

```text
lamin evolve destroy [OPTIONS] TARGETS...

Options:
  -i, --instance TEXT  LaminDB instance to act against, as 'account/instance'.
                       [required]
  --permanent          Delete irreversibly instead of moving objects to the
                       trash.
  --force              Permanently delete the blocking references first, then
                       the targets. This destroys data that cannot be
                       recovered.
  --trash-references   Instead of refusing, move the blocking references that
                       support it to the trash (recoverable) and trash the
                       targets too. Link rows that cannot be trashed are left
                       in place. Cannot be combined with --force or
                       --permanent.
  -y, --yes            Skip the confirmation prompt.
  -c, --config PATH    Path to a lamin_schemas.yml: the destroy's run is
                       recorded under the instance's
                       'track_space'/'track_prefix' there, as 'apply -c' does
                       (falling back to the LAMIN_SCHEMA_EVOLUTION_TRACK_SPACE
                       /LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env vars).
  --help               Show this message and exit.

Arguments:
  TARGETS  One or more plans to act on, each as
           'path/to/plan.py:plan_variable' (the ':plan_variable' suffix is
           optional when the file defines exactly one plan). When several are
           given they are evaluated together: an object dropped from one plan
           but still declared by another listed plan is treated as moved
           between plans, not deleted.
```

::::

::::{dropdown} lamin evolve introspect
:name: cli-evolve-introspect

Generate a plan file from an existing instance.

```text
lamin evolve introspect [OPTIONS] NAME

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  -d, --output-dir PATH     Directory to write the plan file into (default:
                            the current directory).
  -s, --space TEXT          Restrict introspection to one or more space uids
                            (repeat to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from introspection
                            (repeat to add spaces); applied after --space.
  --overwrite               Overwrite the output file if it already exists.
  --no-verify               Skip verifying that applying the plan would make
                            no changes.
  --help                    Show this message and exit.

Arguments:
  NAME  Plan name; the file is written as 'plan_<name>.py'.
```

::::

::::{dropdown} lamin evolve repair
:name: cli-evolve-repair

Repair an instance: dangling objects, duplicates, broken dtypes, stranded values, sheet schemas and misaligned uids.

```text
lamin evolve repair [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  align-uids       Rewrite plan-managed uids on other instances to the...
  dangling         Trash live objects that still reference deleted...
  dtypes           Repair features whose stored dtype string is broken.
  duplicates       Trash duplicate features and schemas, keeping the most...
  sheet-schemas    Detach record types whose schema is a schema type.
  stranded-values  Delete values left behind by a feature that was deleted.
```

::::

::::{dropdown} lamin evolve repair dangling
:name: cli-evolve-repair-dangling

Trash live objects that still reference deleted (trashed) objects.

```text
lamin evolve repair dangling [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  --permanent               Delete irreversibly instead of moving objects to
                            the trash.
  -s, --space TEXT          Restrict pruning to one or more space uids (repeat
                            to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from pruning
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve repair duplicates
:name: cli-evolve-repair-duplicates

Trash duplicate features and schemas, keeping the most recent of each.

```text
lamin evolve repair duplicates [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  --permanent               Delete irreversibly instead of moving objects to
                            the trash.
  -s, --space TEXT          Restrict the scan to one or more space uids
                            (repeat to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from the scan
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve repair dtypes
:name: cli-evolve-repair-dtypes

Repair features whose stored dtype string is broken.

```text
lamin evolve repair dtypes [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  -s, --space TEXT          Restrict verification to one or more space uids
                            (repeat to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from verification
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve repair stranded-values
:name: cli-evolve-repair-stranded-values

Delete values left behind by a feature that was deleted.

```text
lamin evolve repair stranded-values [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  -s, --space TEXT          Restrict the scan to one or more space uids
                            (repeat to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from the scan
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve repair sheet-schemas
:name: cli-evolve-repair-sheet-schemas

Detach record types whose schema is a schema type.

```text
lamin evolve repair sheet-schemas [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  -s, --space TEXT          Restrict the scan to one or more space uids
                            (repeat to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from the scan
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve repair align-uids
:name: cli-evolve-repair-align-uids

Rewrite plan-managed uids on other instances to the source instance's.

```text
lamin evolve repair align-uids [OPTIONS] [TARGETS]...

Options:
  --source TEXT         Instance whose uids win, as 'account/instance'. Its
                        lock files say which uids those are, and are checked
                        against it before anything changes. It is only read.
                        [required]
  -i, --instance TEXT   Target instance to rewrite (repeat for several). With
                        -c and no -i, every configured instance sharing plans
                        with the source.
  -c, --config PATH     lamin_schemas.yml to select targets and plans from,
                        and to read each target's 'track_space' from.
  --merge-duplicates    When a live record already holds the uid an object
                        needs and is the same object (same kind, name, parent
                        and shape -- e.g. a stub 'lamin io sync' created),
                        merge it into the plan's object: move everything
                        referencing it over, then delete it permanently.
                        Likewise, when a schema no plan manages already holds
                        the content of a plan's schema and is not an inferred
                        schema (which is absorbed anyway), merge it into the
                        plan's schema: move everything referencing it over,
                        then move it to the trash.
  --reassign-colliding  When a live record that cannot be merged holds the uid
                        an object needs, give it a fresh uid instead; its
                        references follow it.
  -y, --yes             Skip the per-instance confirmation prompt.
  --dry-run             Stop after each target's rolled-back preview: never
                        prompt, change no instance or lock file and record no
                        run. Lock files that are not committed are only noted.
                        Exits 0 even when a real run would refuse. Cannot be
                        combined with --yes.
  --help                Show this message and exit.

Arguments:
  TARGETS  Plans to align, as 'path/to/plan.py[:plan_variable]'. Their
           prerequisite plans are always included. Omit with -c to align every
           plan the source shares with each target.
```

::::

::::{dropdown} lamin evolve prune-dangling
:name: cli-evolve-prune-dangling

Deprecated alias for 'repair dangling'.

```text
lamin evolve prune-dangling [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  --permanent               Delete irreversibly instead of moving objects to
                            the trash.
  -s, --space TEXT          Restrict pruning to one or more space uids (repeat
                            to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from pruning
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve remove-duplicates
:name: cli-evolve-remove-duplicates

Deprecated alias for 'repair duplicates'.

```text
lamin evolve remove-duplicates [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  --permanent               Delete irreversibly instead of moving objects to
                            the trash.
  -s, --space TEXT          Restrict the scan to one or more space uids
                            (repeat to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from the scan
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve verify-dtypes
:name: cli-evolve-verify-dtypes

Deprecated alias for 'repair dtypes'.

```text
lamin evolve verify-dtypes [OPTIONS]

Options:
  -i, --instance TEXT       LaminDB instance to act against, as
                            'account/instance'.  [required]
  -s, --space TEXT          Restrict verification to one or more space uids
                            (repeat to add spaces).
  -x, --exclude-space TEXT  Exclude one or more space uids from verification
                            (repeat to add spaces); applied after --space.
  -y, --yes                 Skip the confirmation prompt.
  --dry-run                 Print the report and stop: never prompt, change
                            nothing and record no run. Exits 0 even when a
                            real run would refuse. Cannot be combined with
                            --yes.
  -c, --config PATH         Path to a lamin_schemas.yml: the repair's run is
                            recorded under the instance's
                            'track_space'/'track_prefix' there, as 'apply -c'
                            does (falling back to the LAMIN_SCHEMA_EVOLUTION_T
                            RACK_SPACE/LAMIN_SCHEMA_EVOLUTION_TRACK_PREFIX env
                            vars).
  --help                    Show this message and exit.
```

::::

::::{dropdown} lamin evolve utils
:name: cli-evolve-utils

Miscellaneous helpers for authoring plans.

```text
lamin evolve utils [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  add-uid            Add a transform stem uid to an existing plan file.
  fingerprint        Fingerprint one or more plans to prove a refactor is...
  generate-stem-uid  Generate one or more fresh 12-character base62...
```

::::

::::{dropdown} lamin evolve utils add-uid
:name: cli-evolve-utils-add-uid

Add a transform stem uid to an existing plan file.

```text
lamin evolve utils add-uid [OPTIONS] PLAN

Options:
  -i, --instance TEXT  Instance(s) to check the generated uid against for
                       collisions. Repeat to require the uid to be free on
                       several instances.
  --help               Show this message and exit.

Arguments:
  PLAN  Existing plan file to add a PLAN_TRANSFORM_STEM_UID to.
```

::::

::::{dropdown} lamin evolve utils generate-stem-uid
:name: cli-evolve-utils-generate-stem-uid

Generate one or more fresh 12-character base62 transform stem uids.

```text
lamin evolve utils generate-stem-uid [OPTIONS]

Options:
  -i, --instance TEXT        Instance(s) to check the generated uid(s) against
                             for collisions. Repeat to require the uids to be
                             free on several instances.
  -n, --count INTEGER RANGE  Number of stem uids to generate.  [default: 1;
                             x>=1]
  --help                     Show this message and exit.
```

::::

::::{dropdown} lamin evolve utils fingerprint
:name: cli-evolve-utils-fingerprint

Fingerprint one or more plans to prove a refactor is drift-free.

```text
lamin evolve utils fingerprint [OPTIONS] TARGETS...

Options:
  -b, --baseline PATH  Compare the plans against this captured fingerprint and
                       exit non-zero on drift, instead of capturing a new one.
  -o, --output PATH    Write the captured fingerprint JSON here (printed to
                       stdout when omitted); also honored in --baseline mode
                       to save the current one.
  --help               Show this message and exit.

Arguments:
  TARGETS  One or more plans to fingerprint, each as
           'path/to/plan.py:plan_variable' (the ':plan_variable' suffix is
           optional when the file defines exactly one plan). Pass every plan
           the refactor involves: their fingerprints are combined into one.
```

::::

::::{dropdown} lamin evolve ci
:name: cli-evolve-ci

Run plan checks in CI from a lamin_schemas.yml configuration.

```text
lamin evolve ci [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  apply            Apply the configured plans to their instances from the...
  drift            Detect drift across the configured instances for...
  missing-exports  Report the full closure of tracked objects that plans...
  plan             Validate the configured plans and diff them against...
  verify-config    Verify every committed plan is referenced or...
```

::::

::::{dropdown} lamin evolve ci plan
:name: cli-evolve-ci-plan

Validate the configured plans and diff them against their instances.

```text
lamin evolve ci plan [OPTIONS]

Options:
  -c, --config PATH               Path to the CI configuration file (default:
                                  lamin_schemas.yml).  [default:
                                  lamin_schemas.yml]
  --changed TEXT                  Restrict the check to these changed plan
                                  files (repeat to add files).
  -i, --instance TEXT             Restrict the run to this configured instance
                                  (repeat for several); intersects with
                                  --changed.
  -o, --output PATH               Write the Markdown report to this file
                                  (always also printed).
  --fail-on-diff                  Exit non-zero when any instance would see
                                  changes, not only on errors.
  --strict-uids / --no-strict-uids
                                  Count an object whose uid differs from the
                                  uid the plan's lock files for other
                                  instances record as pending, since apply
                                  refuses it (default). --no-strict-uids, or
                                  'strict_uids: false' in the config, only
                                  reports it, matching what 'ci apply' then
                                  does.  [default: strict-uids]
  --strict-renames / --no-strict-renames
                                  Refuse to apply a delete+create that looks
                                  like a rename not declared with ctx.moved
                                  (default). Use --no-strict-renames to allow
                                  it when the objects are genuinely distinct.
                                  'ci plan' reports it as an error for that
                                  plan and instance, failing the run as 'ci
                                  apply' would fail; --no-strict-renames, or
                                  'strict_renames: false' in the config, only
                                  warns.  [default: strict-renames]
  --allow-data-validated-changes  Allow a blocked validation-affecting change
                                  to a feature when the live data proves it
                                  safe: the feature has no stored values, or
                                  every value validates against the new dtype.
                                  Off by default.
                                  'allow_data_validated_changes: true' in the
                                  config opts in too; either way the diff is
                                  the one 'ci apply' applies with the same
                                  setting.
  --preview-migrations            Rehearse each instance's configured plans,
                                  dependencies-first as 'ci apply' runs them,
                                  in one transaction that is always rolled
                                  back: run every declared ctx.migrate()
                                  migration and its checks, and report the
                                  outcome under each diff -- or record the
                                  failure for that plan and instance, as 'ci
                                  apply' would fail. The rehearsal attempts
                                  writes, so the CI account needs write access
                                  to the instances; nothing is persisted and
                                  no lock file is written.
                                  'preview_migrations: true' in the config
                                  opts in too.
  --help                          Show this message and exit.
```

::::

::::{dropdown} lamin evolve ci apply
:name: cli-evolve-ci-apply

Apply the configured plans to their instances from the configuration.

```text
lamin evolve ci apply [OPTIONS]

Options:
  -c, --config PATH               Path to the CI configuration file (default:
                                  lamin_schemas.yml).  [default:
                                  lamin_schemas.yml]
  --changed TEXT                  Restrict the apply to the instances that
                                  have at least one of these changed plan
                                  files (repeat to add files).
  -i, --instance TEXT             Restrict the run to this configured instance
                                  (repeat for several); intersects with
                                  --changed.
  --strict-renames / --no-strict-renames
                                  Refuse to apply a delete+create that looks
                                  like a rename not declared with ctx.moved
                                  (default). Use --no-strict-renames to allow
                                  it when the objects are genuinely distinct.
                                  'strict_renames: false' in the config allows
                                  it too.  [default: strict-renames]
  --allow-data-validated-changes  Allow a blocked validation-affecting change
                                  to a feature when the live data proves it
                                  safe: the feature has no stored values, or
                                  every value validates against the new dtype.
                                  Off by default.
                                  'allow_data_validated_changes: true' in the
                                  config opts in too.
  --strict-uids / --no-strict-uids
                                  Refuse to apply while an object's uid
                                  differs from the uid the plan's lock files
                                  for other instances record (default). --no-
                                  strict-uids, or 'strict_uids: false' in the
                                  config, applies anyway.  [default: strict-
                                  uids]
  -f, --format [terminal|markdown|json]
                                  Output format: 'terminal' for colored output
                                  (default), 'markdown' for a GitHub-flavored
                                  fragment to post on a pull request, or
                                  'json' for a machine-readable document to
                                  pipe into a script.  [default: terminal]
  -o, --output PATH               Also write the rendered document to this
                                  file. Needs --format markdown or --format
                                  json.
  --help                          Show this message and exit.
```

::::

::::{dropdown} lamin evolve ci drift
:name: cli-evolve-ci-drift

Detect drift across the configured instances for scheduled runs.

```text
lamin evolve ci drift [OPTIONS]

Options:
  -c, --config PATH           Path to the CI configuration file (default:
                              lamin_schemas.yml).  [default:
                              lamin_schemas.yml]
  -i, --instance TEXT         Restrict the run to this configured instance
                              (repeat for several).
  -d, --output-dir PATH       Directory to write each instance's drift plan
                              into (default: 'drift_output_dir' from the
                              config, else the current dir).
  -o, --output PATH           Write the Markdown drift report to this file
                              (always also printed).
  --exit-zero                 Always exit 0, even when drift is detected.
  --format-max-chars INTEGER  Truncate the rendered report to at most this
                              many characters (e.g. to fit GitHub's issue-body
                              limit).
  --github-issue              Deprecated: use '--format-max-chars'. Condense
                              the report to fit GitHub's issue-body character
                              limit.
  --help                      Show this message and exit.
```

::::

::::{dropdown} lamin evolve ci missing-exports
:name: cli-evolve-ci-missing-exports

Report the full closure of tracked objects that plans must export.

```text
lamin evolve ci missing-exports [OPTIONS]

Options:
  -c, --config PATH             Path to the CI configuration file (default:
                                lamin_schemas.yml).  [default:
                                lamin_schemas.yml]
  --changed TEXT                Restrict to the instances that declare these
                                plan files (repeat).
  -i, --instance TEXT           Restrict the run to this configured instance
                                (repeat for several); intersects with
                                --changed.
  -o, --output PATH             Write the report to this file (always also
                                printed).
  -f, --format [markdown|json]  Output format: 'markdown' (default) for a
                                report, or 'json' payload.  [default:
                                markdown]
  --json                        Deprecated: use '--format json'. Emit a JSON
                                payload.
  --fail-on-missing             Exit non-zero when any object still needs
                                exporting.
  --help                        Show this message and exit.
```

::::

::::{dropdown} lamin evolve ci verify-config
:name: cli-evolve-ci-verify-config

Verify every committed plan is referenced or explicitly ignored.

```text
lamin evolve ci verify-config [OPTIONS]

Options:
  -c, --config PATH  Path to the CI configuration file (default:
                     lamin_schemas.yml).  [default: lamin_schemas.yml]
  -o, --output PATH  Write the Markdown report to this file (always also
                     printed).
  --help             Show this message and exit.
```

::::

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

The `lamindb` [skill](https://github.com/laminlabs/lamindb/tree/main/lamindb/.agents) ships with the package. After installing `lamindb`, run `uvx library-skills --skill lamindb` so your agent can read it. It will call:

```
lamin track claude   # or: lamin track copilot, lamin track cursor, or lamin track codex
# work with the agent
lamin finish
```

The report includes thinking when the agent stored it as readable text; encoded thoughts are omitted.

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
  agentcore  Start or resume tracking an AgentCore chat.
  claude     Start tracking a Claude Code session in LaminDB.
  codex      Start or resume tracking a Codex session in LaminDB.
  copilot    Start tracking a GitHub Copilot session in LaminDB.
  cursor     Start or resume tracking a Cursor IDE Agent session in LaminDB.
```

→ Python/R alternative: {func}`~lamindb.track` and {func}`~lamindb.finish` for (non-shell) scripts or notebooks

### finish

Finish a tracked session.

This can be a shell script run, or a Claude Code, Copilot, Cursor, or Codex session.

Options:

```text
lamin finish [OPTIONS]

Options:
  --help  Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.finish` for (non-shell) scripts or notebooks

## Administer

### settings

Manage development, cache, modules, branch, space, mount, and run settings.

Get or set the following settings:

- `dev-dir` → development directory {attr}`~lamindb.setup.core.SetupSettings.dev_dir`
- `cache-dir` → cache directory {attr}`~lamindb.setup.core.SetupSettings.cache_dir`
- `modules` → environment schema modules {attr}`~lamindb.setup.core.SetupSettings.modules`
- `branch` → current {attr}`~lamindb.setup.core.SetupSettings.branch`
- `space` → current {attr}`~lamindb.setup.core.SetupSettings.space`
- `mount` → read-only mounts of storage locations, used by `lamin run` to read inputs in place
- `run-where` → default place for `lamin run` (overridden by `$LAMIN_RUN_WHERE` and `--where`)
- `run-access` → default access methods for inputs of `lamin run` (overridden by `$LAMIN_RUN_ACCESS` and `--access`)

You can display your current settings by running: `lamin info`

Examples:

```
# dev-dir (created by lamin init or lamin connect <account/name> --here)
lamin settings dev-dir get
lamin settings dev-dir find .
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
# mount
lamin settings mount storage ./mnt
lamin settings mount path lamin://acme/data/artifact/key/my_file.parquet
lamin settings mount unset ./mnt
# run-where
lamin settings run-where get
lamin settings run-where set modal
lamin settings run-where unset
# run-access
lamin settings run-access get
lamin settings run-access set mount,presigned,credentials
lamin settings run-access unset
```

Options:

```text
lamin settings [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  cache-dir   Get, set, reset, or clear the cache directory.
  dev-dir     Get or find development directories.
  modules     Get or set environment schema modules.
  mount       Mount storage locations read-only via an installed FUSE...
  run-access  Get or set how `lamin run` hands inputs to the target by...
  run-where   Get or set where `lamin run` runs by default.
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

Transfer, import, and export data.

Use `lamin io sync` to sync objects to the current database:

```
lamin io sync https://lamin.ai/laminlabs/lamindata/record/UrcIKR8v0ywim0pE
lamin io sync https://lamin.ai/laminlabs/lamindata/record/UrcIKR8v0ywim0pE --depth 1
lamin io sync https://lamin.ai/laminlabs/lamindata/record/UrcIKR8v0ywim0pE --depth 1 --limit none
lamin io sync record --uid UrcIKR8v0ywim0pE --from laminlabs/lamindata
lamin io sync artifact --key example_datasets/mini_immuno/dataset1.h5ad --from laminlabs/lamindata
```

→ Guide: {doc}`transfer`

Options:

```text
lamin io [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  exportdb  Export registries to parquet files.
  importdb  Import registries from parquet files.
  snapshot  Create an SQLite snapshot of the current database.
  sync      Sync an object to the current database.
```

→ Python/R alternative: {func}`~lamindb.core.sync`

Use `lamin io snapshot` to create an SQLite snapshot of the current database:

```
lamin io snapshot
```

Use `lamin io importdb` to import data from parquet files:

```
lamin io importdb --input-dir path/to/parquet/files
```

Use `lamin io exportdb` to export data to parquet files:

```
lamin io exportdb --output-dir path/to/parquet/files
```

### integrations

Use integrations.

Examples:

```
lamin integrations notion sync db7c1d2ec3a6495e859f8d21d533dd27
lamin integrations notion sync db7c1d2ec3a6495e859f8d21d533dd27 --depth 1 --space my-space --apply
```

Options:

```text
lamin integrations [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  notion  Sync from Notion.
```

→ Python/R alternative: {func}`~lamindb.integrations.notion.sync_objects_from_notion`

## Authenticate

### login

Log into the hub.

`lamin login` prompts for your API key unless you set it via environment variable `LAMIN_API_KEY`.

You can create your API key in your hub account settings at [lamin.ai/settings](https://lamin.ai/settings).

After authenticating once, you can re-authenticate and switch between accounts via `lamin login myhandle`.

Options:

```text
lamin login [OPTIONS] [USER]

Options:
  --help  Show this message and exit.
```

→ Python/R alternative: {func}`~lamindb.setup.login`

### logout

Log out of the hub.

Options:

```text
lamin logout [OPTIONS]

Options:
  --help  Show this message and exit.
```

## Experimental

### run

Run a script or executable, tracked as a run.

Arguments for the target go after `--`. Any `lamin://` URI among them, or the
target itself, is replaced by a local path: read in place from a mounted storage
location if there is one (see `lamin settings mount`), otherwise from the cache.

```
lamin run train.py -- --data lamin://acme/data/artifact/key/train.parquet --epochs 3
lamin run samtools -- view -b lamin://acme/data/artifact/3TrLu3AbQx9dZq2K -o out.bam
lamin run --register-output out.bam align.sh -- --out out.bam
lamin run --where modal --project my_project my_script.py
lamin run --dry-run align.sh -- --out out.bam
```

`--access` lists, in order, which methods may hand an input to the target. Each
URI gets the first one that can serve it, and the log shows which one did and why
earlier ones could not:

- `mount`: a path inside a mount of its storage location, or in local storage
- `cache`: a path in the cache, downloading the artifact if needed
- `presigned`: a time-limited https URL (see `--presign-expiry`); files only
- `credentials`: the unsigned storage URL, e.g. `s3://bucket/key` or
  `gs://bucket/key`, with credentials for the target: `AWS_*` variables (s3) or
  an access token in `CLOUDSDK_AUTH_ACCESS_TOKEN` and `GOOGLE_OAUTH_ACCESS_TOKEN`
  (gs)

The default is `mount,cache`, so a target always gets local paths unless URLs are
allowed. If no allowed method can serve a URI, nothing runs and the error lists
why each failed. The list is as much about what to forbid as about order:
`mount,presigned` never downloads, `mount,cache` never hands out URLs. Set a
default with `$LAMIN_RUN_ACCESS` or `lamin settings run-access`.

A process holds one AWS identity and one GCS token, so `credentials` serves
inputs in argument order until one needs other credentials than an earlier one;
that input falls through to the next allowed method. Mounts and presigned URLs
carry per-location authorization, so allow one of them when inputs live in
storage locations with different credentials.

Environment variables are fixed when the target starts. For long runs on s3, pass
`--credentials-via file`: lamin writes an AWS credentials file, refreshes it from
LaminHub while the target runs, and points `AWS_SHARED_CREDENTIALS_FILE` at it,
which helps tools that re-read it. AWS SDKs (boto3, aws cli) cache what they
read, so use `--credentials-via process` for them instead: `AWS_CONFIG_FILE` then
names a profile whose `credential_process` runs lamin to fetch fresh credentials
whenever the SDK needs them. Tools with their own S3 client that read neither
need `env`, or a mount (`lamin settings mount`).

```
lamin run --access mount,presigned samtools -- view lamin://acme/data/artifact/key/a.bam
lamin run --access credentials aws -- s3 cp lamin://acme/data/artifact/3TrLu3AbQx9dZq2K -
lamin run --access mount,presigned,credentials,cache tool.sh -- lamin://... lamin://...
```

URIs take two forms. The uid form matches nf-lamin; the key form accepts
`space`, `branch` and `version` qualifiers:

```
lamin://<owner>/<instance>/artifact/<uid>[/<subpath>]
lamin://<owner>/<instance>/artifact/key/<key>[/<subpath>][?space=&branch=&version=]
```

Resolved artifacts are linked as run inputs. The target receives
`LAMIN_INITIATED_BY_RUN_UID`, so a script calling `ln.track()` records its run as
a child of this one, plus `LAMIN_INPUT_PATHS` and `LAMIN_MOUNTS` (JSON) to
translate URIs it reads from elsewhere, e.g. config files.

The target's `--version` output, when it has one, is recorded in
`run.params["tool_version"]`, not on its (possibly reused) transform. For a
Python script, the interpreter's `pip freeze` is additionally snapshotted and
linked as `run.environment`, mirroring what `ln.track()` does for its own
process.

Unlike inputs, which carry their own `?branch=&space=` in the URI, the
transform/run/outputs otherwise inherit whatever branch/space the local machine
is pointed at; pass `--branch`/`--space` to make this explicit instead.

An output only registers if the target exits 0, and is uploaded or kept local
following the instance's `keep_artifacts_local` setting; pass `--upload-outputs`
to force it for this run. Use `--dry-run` to preview all of this first.

Options:

```text
lamin run [OPTIONS] TARGET

Options:
  --where [local|modal]           Where to run. Defaults to $LAMIN_RUN_WHERE,
                                  then `lamin settings run-where`, then local.
  --project TEXT                  A valid project name or uid to link the run
                                  to. On Modal, also names the app.
  --register-output TEXT          Register a file the target writes as an
                                  output artifact. Repeatable.
  --access TEXT                   Ordered, comma-separated access methods for
                                  `lamin://` URIs in the arguments; each URI
                                  gets the first that can serve it: `mount`
                                  (path in a mount or local storage), `cache`
                                  (path in the cache, downloading),
                                  `presigned` (https URL), `credentials`
                                  (storage URL plus credentials). Defaults to
                                  $LAMIN_RUN_ACCESS, then `lamin settings run-
                                  access`, then mount,cache.
  --presign-expiry INTEGER RANGE  Seconds a presigned URL stays valid. Only
                                  when --access allows presigned.  [default:
                                  3600; x>=1]
  --credentials-via [env|file|process]
                                  How credentials reach the target when
                                  --access allows credentials: `env`
                                  (environment variables, fixed at launch),
                                  `file` (an AWS credentials file that lamin
                                  keeps fresh), `process` (an AWS config whose
                                  credential_process fetches fresh
                                  credentials; AWS SDKs only). Use file or
                                  process for long runs. s3 only for file and
                                  process.  [default: env]
  --remount                       Remount a storage location if it serves
                                  stale metadata for an input.
  --image-url TEXT                Modal only: a URL to the base docker image.
  --packages TEXT                 Modal only: a comma-separated list of
                                  additional packages.
  --cpu FLOAT                     Modal only: CPU configuration.
  --gpu TEXT                      Modal only: the type of GPU (cuda images
                                  only).
  --branch TEXT                   A branch name or uid for the transform, run
                                  and any outputs. Defaults to the current
                                  local branch.
  --space TEXT                    A space name or uid for the transform, run
                                  and any outputs. Defaults to the current
                                  local space.
  --upload-outputs                Force outputs to upload even if the instance
                                  keeps artifacts local by default.
  --dry-run                       Report what would be linked as inputs and
                                  registered as outputs, without executing the
                                  target or saving anything. Local only.
  --help                          Show this message and exit.
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
