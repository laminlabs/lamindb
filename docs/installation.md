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
jupyter   # parse Jupyter notebooks
bionty    # manage basic biological entities
# cloud backends
aws       # AWS (s3fs, etc.)
gcp       # Google Cloud (gcfs, etc.)
# biological file formats
fcs       # manage FCS files (flow cytometry)
# storage backends
zarr      # store & stream arrays with zarr
# database backends
postgres  # postgres tools
# others
erdiagram # visualize ER diagrams, also needs graphviz
```
