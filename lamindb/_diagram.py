import biogram
from sqlalchemy import MetaData

from lamindb._db._sqlite import get_engine


def diagram(view=True):
    """Make a diagram of entity relationships."""
    engine = get_engine()
    metadata = MetaData(bind=engine)
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
