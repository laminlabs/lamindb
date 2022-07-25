import biogram
import sqlalchemy as sql
from lndb_cli import load_or_create_instance_settings


class schema:
    """Inspect the schema."""

    @classmethod
    def diagram(cls, view=True):
        """Make a diagram of entity relationships."""
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()
        metadata = sql.MetaData(bind=engine)
        metadata.reflect()
        graph = biogram.create_schema_graph(
            metadata=metadata,
            show_datatypes=False,
            show_indexes=False,  # ditto for indexes
            rankdir="TB",
            concentrate=False,  # Don't try to join the relation lines together
        )
        if view:
            biogram.view(graph)
        else:
            return graph

    @classmethod
    @property
    def entities(cls):
        """Return all entities in the db."""
        metadata = sql.MetaData()
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()
        metadata.reflect(bind=engine)
        table_names = [table.name for table in metadata.sorted_tables]
        return table_names

    @classmethod
    @property
    def changes(cls):
        """Return schema changes."""
        return None
