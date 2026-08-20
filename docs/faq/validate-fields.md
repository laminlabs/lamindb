---
execute_via: python
---

# Django field validation

[Django field validation](https://docs.djangoproject.com/en/5.1/ref/validators/) are enabled for models that inherit the `ValidateFields` class.

```bash
# pip install lamindb
mkdir test-django-validation && cd test-django-validation && lamin init
```

```python
import lamindb as ln
from lamindb.core.exceptions import FieldValidationError
```

```python
try:
    ln.Reference(name="my ref", doi="abc.ef", url="myurl.com")
except FieldValidationError as e:
    print(e)
```

```bash
lamin delete --force test-django-validation
```
