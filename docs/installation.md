![pyversions](https://img.shields.io/pypi/pyversions/lamindb)

```shell
pip install lamindb  # basic data management
```

You can configure the installation using `extras`, e.g.,

```shell
pip install 'lamindb[jupyter,bionty]'
```

Supported `extras` are:

```yaml
# commonly used
jupyter   # parse Jupyter notebook metadata
bionty    # basic biological ontologies
# cloud backends
aws       # AWS (s3fs, etc.)
gcp       # Google Cloud (gcfs, etc.)
# biological artifact formats
fcs       # FCS artifacts (flow cytometry)
# storage backends
zarr      # store & stream arrays with zarr
# others
erdiagram # display schema graphs
```

If you'd like a docker container, here is a way: [github.com/laminlabs/lamindb-docker](https://github.com/laminlabs/lamindb-docker).
