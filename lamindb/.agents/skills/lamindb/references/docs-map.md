# LaminDB guide

This skill is the procedure for tracking a session and curating data. The library reference is the guide installed with the package. Open it from the installed package, which does not require a network connection:

```python
from pathlib import Path
import lamindb

print(Path(lamindb.__file__).resolve().parent / ".agents" / "docs")
```

In a git checkout of lamindb, the same pages are in `docs/` at the repository root. Read that tree when you are editing the repository.

Start with these pages:

- `guide.md`
- `tutorial.md`
- `setup.md`
- `query-search.md`
- `track.md`
- `curate.md`
- `transfer.md`

API pages are not in that copy. Use `help()` and the installed Python source.

Use-case, Hub, and pipeline docs are not in the wheel. When the network is available, their index is https://docs.lamin.ai/llms.txt.
