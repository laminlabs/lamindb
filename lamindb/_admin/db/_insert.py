import sqlalchemy as sql

from . import get_engine
from .id import id_file, id_user  # noqa


class insert_if_not_exists:
    """Insert data if it does not yet exist."""

    @classmethod
    def user(cls, user_name):
        from lamindb import load

        df_user = load("user")

        if user_name in df_user.name.values:
            user_id = df_user.index[df_user.name == user_name][0]
            print(f"user {user_name} ({user_id}) already exists")
        else:
            user_id = insert.user(user_name)  # type: ignore
            print(f"added user {user_name} ({user_id})")

        return user_id


class insert:
    """Insert data."""

    @classmethod
    def user(cls, user_name):
        """User."""
        engine = get_engine()
        metadata = sql.MetaData()

        user = sql.Table(
            "user",
            metadata,
            autoload_with=engine,
        )

        with engine.begin() as conn:
            stmt = sql.insert(user).values(name=user_name)
            result = conn.execute(stmt)

        return result.inserted_primary_key[0]

    @classmethod
    def file(cls, name: str, *, source: str = None, source_name: str = None):
        """Data file with its origin."""
        source_id = source
        engine = get_engine()
        metadata = sql.MetaData()

        source_table = sql.Table(
            "file_source",
            metadata,
            sql.Column("id", sql.String, primary_key=True, default=id_user),
            autoload_with=engine,
        )

        file_table = sql.Table(  # primary key gen does not work with reflecting
            "file",
            metadata,
            sql.Column("id", sql.String, primary_key=True, default=id_file),
            autoload_with=engine,
        )

        if source_id is None:
            from nbproject import meta

            source_id = meta.id
            source_name = meta.title
            source_dependency = meta.dependency
            source_type = "nbproject"

            if source_name is None:
                raise RuntimeError(
                    "Can only ingest from notebook with title. Please set a title!"
                )
        else:
            source_dependency = None
            source_type = "other"

        from lamindb import load
        from lamindb._configuration import user_id, user_name

        df_source = load("file_source")
        if source_id not in df_source.index:
            with engine.begin() as conn:
                stmt = sql.insert(source_table).values(
                    id=source_id,
                    name=source_name,
                    dependency=source_dependency,
                    type=source_type,
                    user=user_id,
                )
                conn.execute(stmt)
                print(
                    f"added source {source_name!r} ({source_id}) by user"
                    f" {user_name} ({user_id})"
                )

        with engine.begin() as conn:
            stmt = sql.insert(file_table).values(
                name=name,
                source=source_id,
            )
            result = conn.execute(stmt)
            file_id = result.inserted_primary_key[0]

        return file_id
