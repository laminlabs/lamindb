import lamindb_setup as setup


def graph():
    """Get diagram of entity relationships as `pydot.Dot` graph object.

    It uses `erdiagram.create_schema_graph`
    """
    import erdiagram

    metadata = get_db_metadata()
    return erdiagram.create_schema_graph(
        metadata=metadata,
        show_datatypes=False,
        show_indexes=False,
        rankdir="TB",
        concentrate=True,
    )


def view():
    """View diagram of entity relationships.

    It displays :func:`~lamindb.schema.graph`.
    """
    import erdiagram

    erdiagram.view(graph())


def get_db_metadata():
    import sqlalchemy as sa

    engine = sa.create_engine(setup.settings.instance.db)
    metadata = sa.MetaData(bind=engine)
    metadata.reflect()
    return metadata
