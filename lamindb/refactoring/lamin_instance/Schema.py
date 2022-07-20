import biogram
import sqlalchemy as sql


class Schema:
    """Inspect the schema."""

    def __init__(self, engine):
        self.engine = engine

    def diagram(self, view=True):
        """Make a diagram of entity relationships."""
        metadata = sql.MetaData(bind=self.engine)
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

    @property
    def entities(self):
        """Return all entities in the db."""
        metadata = sql.MetaData()
        metadata.reflect(bind=self.engine)
        table_names = [table.name for table in metadata.sorted_tables]
        return table_names

    @property
    def changes(self):
        """Return schema changes."""
        return None
