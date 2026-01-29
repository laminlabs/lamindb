---
execute_via: python
---

# Django field validation

[Django field validation](https://docs.djangoproject.com/en/5.1/ref/validators/) are enabled for models that inherit the `ValidateFields` class.

```python
# pip install lamindb
!lamin init --storage ./test-django-validation
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

```python
!lamin delete --force test-django-validation
```
