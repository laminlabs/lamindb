from typing import Dict, Optional, Tuple, Union

import sqlmodel as sqm
from lnschema_core import DObject
from sqlmodel import Session


def dobject_to_sqm(entity: Union[sqm.SQLModel, Tuple[sqm.SQLModel]]):
    def if_dobject(entity):
        if entity.__class__.__name__ == "type" and entity.__name__ == "DObject":
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


def get_session_from_kwargs(kwargs: Dict) -> Optional[Session]:
    # modifies kwargs inplace if they contain a session object
    if "session" in kwargs:
        session = kwargs.pop("session")
        assert isinstance(session, Session)
        return session
    else:
        return None
