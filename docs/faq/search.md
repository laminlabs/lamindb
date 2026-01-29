# How does search work?

```python
from laminci.db import setup_local_test_postgres

pgurl = setup_local_test_postgres()
!lamin init --name benchmark_search --db {pgurl} --modules bionty --storage ./benchmark_search
```

Here we show how to perform text search on `SQLRecord` and evaluate some search queries for the {class}`bionty.CellType` ontology.

```python
import lamindb as ln
import bionty as bt

SEARCH_QUERIES_EXACT = (
    "t cell",
    "stem cell",
    "b cell",
    "regulatory B cell",
    "Be2 cell",
    "adipocyte",
)
SEARCH_QUERIES_CONTAINS = ("t cel", "t-cel", "neural", "kidney", "kidne")
TOP_N = 20

bt.CellType.import_source()
```

```python
ln.Record(name="cat[*_*]").save()
```

## Search the registry

```python
for query in SEARCH_QUERIES_EXACT:
    print("Query:", query)
    qs = bt.CellType.search(query)
    display(qs.to_dataframe())

    assert query.lower() == qs[0].name.lower()
```

```python
for query in SEARCH_QUERIES_CONTAINS:
    print("Query:", query)
    qs = bt.CellType.search(query)
    display(qs.to_dataframe())

    top_record = qs[0]
    query = query.lower()
    assert query in top_record.name.lower() or query in top_record.synonyms.lower()
```

Check escaping of special characters.

```python
assert len(ln.Record.search("cat[")) == 1
```

```python
assert len(ln.Record.search("*_*")) == 1
```

## Search the public ontology

```python
ct_public = bt.CellType.public()

df = ct_public.search("b cell", limit=20)
assert df.iloc[0]["name"] == "B cell"
df
```

```python
!docker stop pgtest && docker rm pgtest
!lamin delete --force benchmark_search
```
