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
bionty    # access basic biological ontologies
# cloud backends
aws       # AWS (s3fs, etc.)
gcp       # Google Cloud (gcfs, etc.)
# biological file formats
fcs       # manage FCS files (flow cytometry)
# storage backends
zarr      # store & stream arrays with zarr
# others
erdiagram # visualize ER diagrams, also needs graphviz
```

If you'd like a docker container, here is a way: [github.com/laminlabs/lamindb-docker](https://github.com/laminlabs/lamindb-docker).
