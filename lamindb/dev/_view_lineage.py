from lnschema_core import ORM


def view_lineage(record: ORM, field: str, distance: int = 100):
    if not hasattr(record, "parents"):
        return NotImplementedError(
            f"Lineage view is not supported for {record.__class__.__name__}!"
        )
    import graphviz

    df_edges = _construct_df_edges(record=record, field=field, distance=distance)

    record_label = record.__getattribute__(field)

    u = graphviz.Digraph(record.id, node_attr={"color": "mediumseagreen"})
    u.node(
        record_label.replace(":", "_"),
        label=record_label,
        style="filled",
        color="mediumseagreen",
    )
    for _, row in df_edges.iterrows():
        u.node(row["source"], label=row["source_label"])
        u.edge(row["source"], row["target"], color="darkslategrey")

    return u


def _get_parents(record: ORM, field: str, distance: int):
    """Recursively get parent records within a distance."""
    model = record.__class__
    condition = f"children__{field}"
    results = model.select(**{condition: record.__getattribute__(field)}).all()
    if distance < 2:
        return results

    for d in range(2, distance):
        condition = "children__" + condition
        records = model.select(**{condition: record.__getattribute__(field)}).all()

        if len(records) == 0:
            return results

        results = results | records
    return results


def _construct_df_edges(record: ORM, field: str, distance: int):
    """Construct a DataFrame of edges as the input of graphviz.Digraph."""
    parents = _get_parents(record=record, field=field, distance=distance)
    records = parents | record.__class__.objects.filter(id=record.id)
    df = records.distinct().df(include=[f"parents__{field}"])
    df_edges = df[[f"parents__{field}", field]]
    df_edges = df_edges.explode(f"parents__{field}")
    df_edges.dropna(axis=0, inplace=True)
    df_edges.rename(
        columns={f"parents__{field}": "source", field: "target"}, inplace=True
    )
    df_edges = df_edges.drop_duplicates()

    # colons messes with the node formatting:
    # https://graphviz.readthedocs.io/en/stable/node_ports.html
    df_edges["source_label"] = df_edges["source"]
    df_edges["target_label"] = df_edges["target"]
    df_edges["source"] = df_edges["source"].str.replace(":", "_")
    df_edges["target"] = df_edges["target"].str.replace(":", "_")
    return df_edges
