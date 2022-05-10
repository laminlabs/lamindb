# lamindb

lamindb is the biological data management system.

## Installation and configuration

Install the development version:
```
flit install -s --deps develop
```
Then configure lamin for use with Notion:
```
lamin configure
```
and AWS (`pip install awscli`):
```
aws configure
```

## Tutorials

* Notion integration: https://github.com/laminlabs/lamin/blob/main/tests/notion.ipynb
* Ingesting Eraslan21: https://github.com/laminlabs/usecases/blob/main/ingest-eraslan21.ipynb
