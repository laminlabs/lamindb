from typing import Tuple, Union

import sqlmodel as sqm
from lndb_setup import settings
from lnschema_core import DObject


def session() -> sqm.Session:
    """Get connection session to DB engine.

    Returns a `sqlmodel.Session` object.
    """
    return settings.instance.session()


def dobject_func_to_class(entity: Union[sqm.SQLModel, Tuple[sqm.SQLModel]]):
    def if_dobject(entity):
        if entity.__class__.__name__ == "function" and entity.__name__ == "DObject":
            return DObject
        else:
            return entity

    if isinstance(entity, tuple):
        entities = list(entity)
        for i, ent in enumerate(entities):
            entities[i] = if_dobject(ent)
        return entities
    else:
        return if_dobject(entity)
