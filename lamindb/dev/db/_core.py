from typing import Dict, Optional, Tuple, Union

import sqlmodel as sqm
from lamindb_setup import _USE_DJANGO
from lnschema_core import File
from sqlmodel import Session


def file_to_sqm(entity: Union[sqm.SQLModel, Tuple[sqm.SQLModel]]):
    def if_file(entity):
        if entity.__class__.__name__ == "type" and entity.__name__ == "File":
            return File
        else:
            return entity

    if isinstance(entity, tuple):
        entities = list(entity)
        for i, ent in enumerate(entities):
            entities[i] = if_file(ent)
        return entities
    else:
        return if_file(entity)


def get_session_from_kwargs(kwargs: Dict) -> Optional[Session]:
    # modifies kwargs inplace if they contain a session object
    if "session" in kwargs:
        session = kwargs.pop("session")

        if _USE_DJANGO:
            from lamindb_setup.dev._settings_instance import DjangoSession

            assert isinstance(session, DjangoSession)
        else:
            from sqlmodel import Session

            assert isinstance(session, Session)
        return session
    else:
        return None
