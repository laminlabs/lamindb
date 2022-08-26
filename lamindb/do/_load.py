import pandas as pd
import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from lnschema_core import id
from nbproject import meta

from ..dev import filepath_from_dobject
from ..dev.file import load_to_memory
from ..schema import core


def populate_dtransform_in(dobject):
    jupynb_id = meta.store.id
    jupynb_v = meta.store.version  # version to be set in publish()
    jupynb_name = meta.live.title
    engine = settings.instance.db_engine()

    with sqm.Session(engine) as session:
        result = session.get(core.jupynb, (jupynb_id, jupynb_v))
        if result is None:
            session.add(
                core.jupynb(
                    id=jupynb_id, v=jupynb_v, name=jupynb_name, user_id=settings.user.id
                )
            )
            dtransform_id = id.id_dtransform()
            session.add(
                core.dtransform(
                    id=dtransform_id,
                    jupynb_id=jupynb_id,
                    jupynb_v=jupynb_v,
                )
            )
            session.commit()
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user.handle}."
            )
        else:
            dtransform_id = session.exec(
                sqm.select(core.dtransform).where(
                    core.dtransform.jupynb_id == jupynb_id,
                    core.dtransform.jupynb_v == jupynb_v,
                )
            ).first()  # change to .one() as soon as dtransform ingestion bug fixed
        session.add(
            core.dtransform_in(
                dtransform_id=dtransform_id,
                dobject_id=dobject.id,
                dobject_v=dobject.v,
            )
        )


class load:
    """Load data."""

    @classmethod
    def entity(cls, entity_name) -> pd.DataFrame:
        """Load observations of entity as dataframe."""
        engine = settings.instance.db_engine()
        with engine.connect() as conn:
            df = pd.read_sql_table(entity_name, conn)
            if "id" in df.columns:
                if "v" in df.columns:
                    df = df.set_index(["id", "v"])
                else:
                    df = df.set_index("id")
        return df

    @classmethod
    def dobject(cls, dobject: core.dobject):
        """Load dobject into memory."""
        filepath = filepath_from_dobject(dobject)
        populate_dtransform_in(dobject)
        return load_to_memory(filepath)
