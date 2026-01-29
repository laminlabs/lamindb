---
executable_via: python
---

# What happens if I import a schema module without lamindb?

```python
# !pip install 'lamindb[bionty]'
!lamin init --storage testmodule --modules bionty
```

Upon `import`, nothing yet happens:

```python
import bionty as bt
```

If you try to access an attribute (other than `model`), you'll load the instance in the same way as calling `import lamindb`.

Under the hood, `lamindb` is imported!

```python
assert bt.Organism(name="human") is not None
```

```python
!lamin delete --force testmodule
```
