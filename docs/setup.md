# Install & setup

To install LaminDB on ![pyversions](https://img.shields.io/pypi/pyversions/lamindb) with data-science-related dependencies, run:

```shell
pip install lamindb
```

For a minimal installation of the `lamindb` namespace, use:

```shell
pip install lamindb-core
```

You can configure the installation using `extras`, e.g.,

```shell
pip install 'lamindb[gcp]'  # install with support for the Google Cloud Platform
```

To see all configuration options, read the `pyproject.toml` [here](https://github.com/laminlabs/lamindb/blob/main/pyproject.toml).
If you'd like to install a development version from GitHub, consider these [notes](https://github.com/laminlabs/lamindb/blob/main/CONTRIBUTING.md).

## Create a database

You can create a LaminDB instance using the [init command](https://docs.lamin.ai/cli#init) with these options:

- `storage`: a default storage location (e.g. `s3://my-bucket`, `gs://my-bucket`, `./my-data-dir`)
- `name` (optional): a name (e.g., `my-assets`)
- `db` (optional): a Postgres database connection URL, do not pass for SQLite
- `modules` (optional): comma-separated string of lamindb modules

If you are only interested in tracking artifacts and their transformations, init your local SQLite database via:

```shell
lamin init --storage ./mydata
```

Mount the Bionty module:

```shell
lamin init --storage mydata --modules bionty
```

You can also pass an AWS S3 bucket:

```shell
lamin init --storage s3://<bucket_name> --modules bionty
```

Instead of SQLite, you can pass a Postgres connection string:

```shell
lamin init --storage gs://<bucket_name> --db postgresql://<user>:<pwd>@<hostname>:<port>/<dbname> --modules bionty
```

To delete a database, call:

```shell
lamin delete instance name
```

This will only work if the database has no data in its storage location.

## Connect to a database

Connect to a database for reads:

```python
import lamindb as ln

db = ln.DB("account/name")
```

Configure your default database on the terminal:

```shell
lamin connect <account/name>  # tip: add flag `--here` to scope to current directory
```

In Python/R, you'll now auto-connect. To disconnect, run `lamin disconnect`.

To configure your default database in a Python/R session, run:

```python
ln.connect("account/name")
```

To access private databases through the hub, you need an account. It's free & [signing up](https://lamin.ai/signup) takes <1 min. To log in, run:

```shell
lamin login
```

You will be prompted for your API key. You can create your API key on your [account settings page](https://lamin.ai/settings).

If you have multiple accounts that already logged into a compute environment, you can switch between them using your handle:

```shell
lamin login mclintock
```

Log out:

```shell
lamin logout
```

## Configure settings

Print info about settings on the terminal:

```shell
lamin info
```

Settings persist in `~/.lamin/` (configurable via `LAMIN_SETTINGS_DIR`) and can also be accessed via Python:

```python
import lamindb as ln

ln.setup.settings
```

This returns a {class}`~lamindb.setup.core.SetupSettings` object.

### Use paths with AWS-S3-compatible endpoints

It is possible to create a database with a path that uses an AWS-S3-compatible endpoint url. Such endpoints allow to access non-S3 buckets using the same API that is used for S3:

```shell
lamin init --storage s3://<bucket_name>?endpoint_url=http://endpoint.com:port
```

This assumes that the endpoint url is `http://endpoint.com` with a port specified.

It is also possible to set a path with s3-compatible endpoint as a default storage for an existing instance for the current python session.

```python
import lamindb as ln

ln.settings.storage = ln.UPath("s3://<bucket_name>", endpoint_url="http://endpoint.com:port")
```

### Manage the cache directory

`lamindb` maintains a local cache for files and folders stored in the cloud (e.g., AWS S3, Google Cloud Storage, HTTP, Hugging Face, etc.).

When an {class}`~lamindb.Artifact` object representing a file or folder in the cloud is accessed for the first time via {meth}`~lamindb.Artifact.cache` or {meth}`~lamindb.Artifact.load`, it is downloaded to the cache. Subsequent accesses read from the cached copy, as long as the original file or folder did not change.

The cache directory can be accessed via {class}`lamindb.settings`:

```python
ln.settings.cache_dir
```

Or via the CLI:

```shell
lamin settings cache-dir get
```

It can be configured via the CLI or by setting the `LAMIN_CACHE_DIR` environment variable. Here is the CLI command:

```shell
lamin settings cache-dir set some/path/to/cache
```

### Configuring a system-wide cache

If you are using `lamindb` on a multi-user system such as a shared compute cluster, you can configure a shared default cache for all users to avoid duplicating cached data for each individual user.

To set this up, first find the location of the `lamindb` system settings directory:

```shell
lamin info
```

In the `Local directories` section, locate the path shown in `system settings` - this is the directory you need. In this directory you need to create a text file `system.env` that contains a line with the path you need for the system cache folder (repalce `absolute/path/to/your/system/cache` with your path):

```
lamindb_cache_path=absolute/path/to/your/system/cache
```

This cache folder will be used by default for all users on the system unless they explicitly configure their own cache folder with CLI `lamin cache set`.

## Database modules

1. Any LaminDB instance can mount custom schema modules with any number of registries
2. Each schema module is a Python package that defines registries using the {class}`~lamindb.models.SQLRecord` class
3. Every registry corresponds to a SQL table in the underlying Postgres or SQLite database

The core database schema is built into the `lamindb` API. Most of LaminDB's central classes ({class}`~lamindb.Artifact`, {class}`~lamindb.Transform`, {class}`~lamindb.User`, etc.) are registries. You can see their source code [here](https://github.com/laminlabs/lnschema-core/blob/main/lamindb/models.py).

### Compatibility matrix

Below is the compatibility matrix for the core `lamindb` schema. To upgrade the state of the SQL database (`DB`) from a lower version to your current Python package (`PP`) version, you call: `lamin migrate deploy`

| --         | PP 2.9 | PP 2.6 | PP 2.4 | PP 2.2 | PP 2.1 | PP 2.0 |
| ---------- | ------ | ------ | ------ | ------ | ------ | ------ |
| **DB 2.9** | 🟢     | 🟢     | 🟢     | 🟢     | 🟢     | 🟢     |
| **DB 2.6** | 🔴     | 🟢     | 🟢     | 🟢     | 🟢     | 🟢     |
| **DB 2.4** | 🔴     | 🔴     | 🟢     | 🟢     | 🟢     | 🟢     |
| **DB 2.2** | 🔴     | 🔴     | 🔴     | 🟢     | 🟢     | 🟢     |
| **DB 2.1** | 🔴     | 🔴     | 🔴     | 🔴     | 🟢     | 🟢     |
| **DB 2.0** | 🔴     | 🔴     | 🔴     | 🔴     | 🟢     | 🟢     |

### Custom modules

You can set up your own modules & registries or [reach out](https://lamin.ai/contact) for support within Lamin's Team or Enterprise plan.

You'll see how simple it is if you look at this example: [pertdb/models.py](https://github.com/laminlabs/pertdb/blob/main/pertdb/models.py). You only need a single Python file to define registries via data models.

If you are an admin, you can use two commands to create and deploy migrations:

- `lamin migrate create` (only needed when creating your own custom modules)
- `lamin migrate deploy`

### Create a migration

You need to have the package installed locally:

```shell
git clone https://github.com/my-org/my-module
cd my-module
pip install -e .
```

Edit the registries in your module.

Then, call

```shell
lamin migrate create
```

to create the migration script.

When you're happy, commit them to your GitHub repo, and ideally make a new release.

### Deploy a migration

To deploy the migration call `lamin migrate deploy`.

<!-- #region -->

## FAQ

### Where is the SQLite file of a LaminDB instance?

The SQLite file is in the default storage location of the instance and called `lamin.db`.

You can also see it as part of the database connection string:

```python
ln.setup.settings.instance.db
#> sqlite:///path-to-sqlite
```

If the default storage is in the cloud, the SQLite file is cached in the local cache directory ({attr}`~lamindb.setup.core.StorageSettings.cache_dir`):

```python
ln.setup.settings.storage.cache_dir
#> path-to-cache-dir
```

<!-- #endregion -->
