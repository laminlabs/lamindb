# `lamin`

## `lamin login`

```
usage: ipykernel_launcher.py login [-h] [--key key] [--password pw] user

positional arguments:
  user           Email or user handle. Email is needed at first login.

options:
  -h, --help     show this help message and exit
  --key key      API key or legacy password.
  --password pw  API key or legacy password.
```

## `lamin init`

```
usage: ipykernel_launcher.py init [-h] [--storage s] [--db d]
                                  [--schema schema] [--name n] [--vault]

options:
  -h, --help       show this help message and exit
  --storage s      Storage root. Either local dir, ``s3://bucket_name`` or
                   ``gs://bucket_name``.
  --db d           Database connection url, do not pass for SQLite.
  --schema schema  Comma-separated string of schema modules. None if not set.
  --name n         Instance name.
  --vault          Use vault to manage credentials.
```

## `lamin load`

```
usage: ipykernel_launcher.py load [-h] [--db d] [--storage s] [--vault] i

positional arguments:
  i            The instance identifier can the instance name (owner is current
               user), handle/name, or the URL: https://lamin.ai/handle/name.

options:
  -h, --help   show this help message and exit
  --db d       Database connection url, do not pass for SQLite.
  --storage s  Load the instance with an updated default storage.
  --vault      Use vault to manage credentials.
```

## `lamin close`

```
usage: ipykernel_launcher.py close [-h]

options:
  -h, --help  show this help message and exit
```

## `lamin delete`

```
usage: ipykernel_launcher.py delete [-h] [--force] i

positional arguments:
  i           Instance name.

options:
  -h, --help  show this help message and exit
  --force     Do not ask for confirmation
```

## `lamin track`

```
usage: ipykernel_launcher.py track [-h] [--pypackage pypackage] filepath

positional arguments:
  filepath              A path to the notebook.

options:
  -h, --help            show this help message and exit
  --pypackage pypackage
                        One or more (delimited by ',') python packages to
                        track.
```

## `lamin info`

```
usage: ipykernel_launcher.py info [-h]

options:
  -h, --help  show this help message and exit
```

## `lamin migrate`

```
usage: ipykernel_launcher.py migrate [-h] [--package-name PACKAGE_NAME]
                                     [--end-number END_NUMBER]
                                     [--start-number START_NUMBER]
                                     {create,deploy,squash}

positional arguments:
  {create,deploy,squash}
                        Manage migrations.

options:
  -h, --help            show this help message and exit
  --package-name PACKAGE_NAME
  --end-number END_NUMBER
  --start-number START_NUMBER
```

## `lamin save`

```
usage: ipykernel_launcher.py save [-h] filepath

positional arguments:
  filepath    A path to the notebook.

options:
  -h, --help  show this help message and exit
```

## `lamin set`

```
usage: ipykernel_launcher.py set [-h] [--storage f]

options:
  -h, --help   show this help message and exit
  --storage f  Storage root. Either local dir, ``s3://bucket_name`` or
               ``gs://bucket_name``.
```

## `lamin schema`

```
usage: ipykernel_launcher.py schema [-h] {view}

positional arguments:
  {view}      View schema.

options:
  -h, --help  show this help message and exit
```
