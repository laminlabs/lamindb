from typing import Union

from django.db.models import Manager, QuerySet
from lnschema_core import BaseORM


def select(*entity: BaseORM, **fields) -> Union[QuerySet, Manager]:
    """Query data.

    Guide: :doc:`/guide/select`.

    Args:
        entity: Table, tables, or tables including column specification.
        fields: Fields and values passed as keyword arguments.

    Returns:
        A QuerySet or Django Manager.
    """
    if len(entity) > 1:
        raise NotImplementedError  # currently no longer implemented with Django
    from lnschema_core.models import LaminQuerySet

    manager = LaminQuerySet.as_manager()
    manager.model = entity[0]
    if len(fields) > 0:
        return manager.filter(**fields)
    else:
        return manager
