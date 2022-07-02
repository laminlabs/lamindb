import lamindb as db
import pandas as pd
import sqlmodel as sqm
from lamindb import setup
from sqlalchemy import create_engine

from ..admin.db import get_engine

import analytics

analytics.write_key = 'aBGwRpa5B1EGbYsQdHOqktau4pFnJdpE'

class EventTacker:

    def __init__(self):
        self.user_db_engine = None
        self.internal_db_engine = None

    def track_ingest(self, file_id):

        from nbproject import meta

        parameters = {
            "type": "ingest",
            "user_id": setup.settings().user_id,
            "interface_id": meta.store.id,
            "file_id": file_id
        }

        user_db_instance = self.__insert_event_into_sql_lite(
            self.__get_user_db_engine(),
            db.model.track_do,
            parameters
        )

        self.__insert_event_into_supabase(
            self.__get_internal_db_engine(),
            "events_ingest",
            parameters
        )

        self.__insert_event_into_segment('Ingest', parameters)

        return user_db_instance

    ### Helpers function to insert data into different kind of databases

    def __insert_event_into_sql_lite(self, engine, sql_lite_model, params):

        with sqm.Session(engine) as session:
            instance = sql_lite_model(**params)
            session.add(instance)
            session.commit()
            session.refresh(instance)

        return instance

    def __insert_event_into_supabase(self, engine, table_name, params):
        df = pd.DataFrame([params])
        df.to_sql(table_name, engine,
                  if_exists="append", index=False)

    def __insert_event_into_segment(self, event_name, parameters):
        analytics.track(setup.settings().user_id, event_name, parameters)

    ### Helpers function to get db engine

    def __get_user_db_engine(self):
        if self.user_db_engine is None:
            self.user_db_engine = get_engine()
        return self.user_db_engine

    def __get_internal_db_engine(self):
        if self.internal_db_engine is None:
            self.internal_db_engine = create_engine(
                "postgresql://postgres:B9N!dajGJ2M2@db.dzmuroeojlpdhyqfwyjt.supabase.co:5432/postgres")
        return self.internal_db_engine
