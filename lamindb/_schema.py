import biogram
import sqlalchemy as sql

from lamindb._admin.db import get_engine


class schema:
    """Inspect the schema."""

    @classmethod
    def diagram(cls, view=True):
        """Make a diagram of entity relationships."""
        engine = get_engine()
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
        engine = get_engine()
        metadata.reflect(bind=engine)
        table_names = [table.name for table in metadata.sorted_tables]
        return table_names

    @classmethod
    @property
    def changes(cls):
        """Return all schema changes."""
        raise NotImplementedError
