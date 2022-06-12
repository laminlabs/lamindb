import sqlalchemy as sql

from ...dev.id import id_file, id_user  # noqa
from . import get_engine


class insert_if_not_exists:
    """Insert data if it does not yet exist."""

    @classmethod
    def user(cls, user_name):
        from lamindb.do import load

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
            sql.Column("id", sql.String, primary_key=True, default=id_user),
            autoload_with=engine,
        )

        with engine.begin() as conn:
            stmt = sql.insert(user).values(name=user_name)
            result = conn.execute(stmt)

        return result.inserted_primary_key[0]

    @classmethod
    def file(cls, name: str, *, interface: str = None, interface_name: str = None):
        """Data file with its origin."""
        interface_id = interface
        engine = get_engine()
        metadata = sql.MetaData()

        interface_table = sql.Table(
            "interface",
            metadata,
            autoload_with=engine,
        )

        file_table = sql.Table(  # primary key gen does not work with reflecting
            "file",
            metadata,
            sql.Column("id", sql.String, primary_key=True, default=id_file),
            autoload_with=engine,
        )

        if interface_id is None:
            from nbproject import meta

            interface_id = meta.id
            interface_name = meta.title
            interface_dependency = meta.dependency
            interface_type = "nbproject"

            if interface_name is None:
                raise RuntimeError(
                    "Can only ingest from notebook with title. Please set a title!"
                )
        else:
            interface_dependency = None
            interface_type = "other"

        from lamindb._configuration import user_id, user_name
        from lamindb.do import load

        df_interface = load("interface")
        if interface_id not in df_interface.index:
            with engine.begin() as conn:
                stmt = sql.insert(interface_table).values(
                    id=interface_id,
                    name=interface_name,
                    dependency=interface_dependency,
                    type=interface_type,
                    user=user_id,
                )
                conn.execute(stmt)
                print(
                    f"added interface {interface_name!r} ({interface_id}) by user"
                    f" {user_name} ({user_id})"
                )

        with engine.begin() as conn:
            stmt = sql.insert(file_table).values(
                name=name,
                interface=interface_id,
            )
            result = conn.execute(stmt)
            file_id = result.inserted_primary_key[0]

        return file_id
