from typing import Union

from django.db.models import Manager
from lnschema_core import ORM

from lamindb._queryset import QuerySet


def select(*ORM: ORM, **expressions) -> Union[QuerySet, Manager]:
    """Query records.

    Guide: :doc:`/guide/select`.

    Args:
        ORM: An ORM class.
        expressions: Fields and values passed as Django query expressions.

    Returns:
        A `QuerySet` or Django `Manager`.
    """
    if len(ORM) > 1:
        raise NotImplementedError(
            "Currently no longer implemented with Django, but we hope to have it back."
        )
    manager = QuerySet.as_manager()
    manager.model = ORM[0]
    if len(expressions) > 0:
        return manager.filter(**expressions)
    else:
        return manager
