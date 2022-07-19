from pathlib import Path
from typing import Union

import pandas as pd
from cloudpathlib import CloudPath
from sqlmodel import Session, SQLModel, create_engine, select


class SQLiteLocalDbClient:
    def __init__(
        self, instance_name: str, db_base_path: Union[Path, CloudPath]
    ) -> None:
        self.instance_name = instance_name
        self.db_base_path = db_base_path.absolute()
        self.db_file_path = db_base_path / f"{instance_name}.db"
        self.engine = create_engine(f"sqlite:///{str(self.db_file_path)}", future=True)

        if not self.db_file_path.exists():
            self.__setup()

    def __setup(self) -> None:
        SQLModel.metadata.create_all(self.engine)

    def insert(self, sqlite_model, params):
        with Session(self.engine) as session:
            instance = sqlite_model(**params)
            session.add(instance)
            session.commit()
            session.refresh(instance)
        return instance

    def load_rows(self, sqlite_model: SQLModel, match_keys):
        with Session(self.engine) as session:
            statement = select(sqlite_model)
            for key, value in match_keys.items():
                statement = statement.where(key == value)
            instances = session.exec(statement).all()
            return instances

    def load_table(self, table_name):
        with self.engine.connect() as conn:
            df = pd.read_sql_table(table_name, conn)
            return df
