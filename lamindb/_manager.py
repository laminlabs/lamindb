from typing import Optional

from django.db import models


class Manager(models.Manager):
    """Extension of Django QuerySet."""

    def list(self, field: Optional[str] = None):
        """Populate a list with the results."""
        if field is None:
            return [item for item in self.all()]
        else:
            return [item for item in self.values_list(field, flat=True)]


setattr(models.Manager, "list", Manager.list)
