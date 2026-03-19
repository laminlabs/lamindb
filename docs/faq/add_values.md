# `ln.Record` features: `add_values` vs `set_values`

---

### `add_values`

- Adds annotations. For categorical features, multiple calls for the same feature name typically append links rather than overwrite.
- For `float` (and similar), there is only one slot per feature per record.

After two `add_values` calls on `cell_line` (`HEK293` then `HeLa cell`), `describe()` shows both values on one feature — accumulation from `add_values`, not a single replacement.

### `set_values`

- Like `add_values`, but first removes all existing feature annotations on that record, then applies the dictionary you pass.
- Use it when you want replacement-style behavior (e.g. only one cell line, or a clean slate).

```python
"""Record features: add_values (append for categorical) vs set_values (replace)."""

import bionty as bt
import lamindb as ln

ln.connect("...")


experiment_1 = ln.Record(name="Experiment1").save()

cell_line = ln.Feature(name="cell_line", dtype=bt.CellLine).save()
gc_content = ln.Feature(name="gc_content", dtype=float).save()
experiment_1.features.add_values({"gc_content": 4.5})
experiment_1.features.add_values({"cell_line": "HEK293"})
experiment_1.features.add_values({"cell_line": "HeLa cell"})

################################################################
Record: Experiment1
├── uid: LLxnTXw9OMML2yLc                run:
│   type: Experiment                     is_type: False
│   schema:                              reference:
│   branch: main                         space: all
│   created_at: 2026-03-19 18:07:59 UTC  created_by: biocoder1001
└── Features
    ├── cell_line         cat[bionty.CellLine]     HEK293, HeLa cell
    └── gc_content        float                    4.5
###############################################################



# experiment_1.features.add_values({"gc_content": 5.6})  # IntegrityError

experiment_1.features.set_values(
    {"cell_line": "HeLa cell"}
)

experiment_1.describe()
###################################################################
Record: Experiment1
├── uid: LLxnTXw9OMML2yLc                run:
│   type: Experiment                     is_type: False
│   schema:                              reference:
│   branch: main                         space: all
│   created_at: 2026-03-19 18:07:59 UTC  created_by: biocoder1001
└── Features
    ├── cell_line         cat[bionty.CellLine]     HeLa cell
    └── gc_content        float                    4.5
###################################################################

---



```
