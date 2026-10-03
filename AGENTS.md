# Docs

Use single backticks for inline code, never double backticks.

## Guide

In this checkout, start at `docs/guide.md` and grep `docs/`.

After `pip install lamindb`, the same guide is at `lamindb/.agents/docs/` inside the installed package:

```python
from pathlib import Path
import lamindb

Path(lamindb.__file__).resolve().parent / ".agents" / "docs"
```

API pages in `docs/` are autodoc stubs. Use `help()` or the Python source. The public index, when the network is available, is https://docs.lamin.ai/llms.txt.
