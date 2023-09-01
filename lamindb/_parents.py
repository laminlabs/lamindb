import builtins
from typing import List, Optional, Set, Union

from lamin_utils import logger
from lnschema_core import Dataset, File, Registry, Run, Transform
from lnschema_core.models import HasParents, format_field_value

from lamindb._utils import attach_func_to_class_method

from . import _TESTING
from ._registry import StrField, get_default_str_field

LAMIN_GREEN_LIGHTER = "#10b981"
LAMIN_GREEN_DARKER = "#065f46"
GREEN_FILL = "honeydew"
TRANSFORM_EMOJIS = {"notebook": "📔", "app": "🖥️", "pipeline": "🧩"}
is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


def _transform_emoji(transform: Transform):
    if transform is not None:
        return TRANSFORM_EMOJIS.get(transform.type, "💫")
    else:
        return TRANSFORM_EMOJIS["pipeline"]


def _view(u):
    try:
        if is_run_from_ipython:
            from IPython.display import display

            return display(u)
        else:
            return u
    except (FileNotFoundError, RuntimeError):  # pragma: no cover
        logger.error(
            "please install the graphviz executable on your system:\n  - Ubuntu: `sudo"
            " apt-get install graphviz`\n  - Windows:"
            " https://graphviz.org/download/#windows\n  - Mac: `brew install graphviz`"
        )


def view_parents(
    self,
    field: Optional[StrField] = None,
    with_children: bool = False,
    distance: int = 5,
):
    if field is None:
        field = get_default_str_field(self)
    if not isinstance(field, str):
        field = field.field.name

    return _view_parents(
        record=self, field=field, with_children=with_children, distance=distance
    )


def view_flow(file: File, with_children: bool = True):
    """Graph of data flow.

    Notes:
        For more info, see use cases: :doc:`docs:data-flow`.

    Examples:
        >>> file.view_flow()
    """
    import graphviz

    all_runs = _get_all_parent_runs(file)
    if with_children:
        all_runs.update(_get_all_child_runs(file))
    df_edges = _df_edges_from_runs(all_runs)

    file_label = _label_file_run_transform(file)

    u = graphviz.Digraph(
        file.id,
        node_attr={
            "fillcolor": GREEN_FILL,
            "color": LAMIN_GREEN_DARKER,
            "fontname": "Helvetica",
            "fontsize": "10",
        },
        edge_attr={"arrowsize": "0.5"},
    )

    def add_node(
        record: Union[Run, File], node_id: str, node_label: str, u: graphviz.Digraph
    ):
        if isinstance(record, Run):
            style = "rounded,filled"
            shape = "box"
            fillcolor = "gainsboro"
        else:
            style = "rounded,filled"
            shape = "box"
            fillcolor = GREEN_FILL
        u.node(
            node_id,
            label=node_label,
            shape=shape,
            style=style,
            fillcolor=fillcolor,
        )

    for _, row in df_edges.iterrows():
        add_node(row["source_record"], row["source"], row["source_label"], u)
        if row["target_record"] not in df_edges["source_record"]:
            add_node(row["target_record"], row["target"], row["target_label"], u)

        u.edge(row["source"], row["target"], color="dimgrey")
    # label the searched file
    u.node(
        file.id,
        label=file_label,
        style="rounded,filled",
        fillcolor=LAMIN_GREEN_LIGHTER,
        shape="box",
    )

    _view(u)


def _view_parents(
    record: Registry, field: str, with_children: bool = False, distance: int = 100
):
    """Graph of parents."""
    if not hasattr(record, "parents"):
        raise NotImplementedError(
            f"Parents view is not supported for {record.__class__.__name__}!"
        )
    import graphviz
    import pandas as pd

    df_edges = _df_edges_from_parents(record=record, field=field, distance=distance)
    if with_children:
        df_edges = pd.concat(
            [
                df_edges,
                _df_edges_from_parents(
                    record=record, field=field, distance=distance, children=True
                ),
            ]
        ).drop_duplicates()

    record_label = record.__getattribute__(field)

    u = graphviz.Digraph(
        record.id,
        node_attr={
            "color": LAMIN_GREEN_DARKER,
            "fillcolor": GREEN_FILL,
            "shape": "box",
            "style": "rounded,filled",
            "fontname": "Helvetica",
            "fontsize": "10",
        },
        edge_attr={"arrowsize": "0.5"},
    )
    u.node(
        record.id,
        label=_label_file_run_transform(record)
        if record.__class__.__name__ == "Transform"
        else _add_emoji(record, record_label),
        fillcolor=LAMIN_GREEN_LIGHTER,
    )
    for _, row in df_edges.iterrows():
        u.node(row["source"], label=row["source_label"])
        u.node(row["target"], label=row["target_label"])
        u.edge(row["source"], row["target"], color="dimgrey")

    _view(u)


def _get_parents(record: Registry, field: str, distance: int, children: bool = False):
    """Recursively get parent records within a distance."""
    if children:
        key = "parents"
    else:
        key = "children"
    model = record.__class__
    condition = f"{key}__{field}"
    results = model.filter(**{condition: record.__getattribute__(field)}).all()
    if distance < 2:
        return results

    d = 2
    while d < distance:
        condition = f"{key}__{condition}"
        records = model.filter(**{condition: record.__getattribute__(field)}).all()

        if len(records) == 0:
            return results

        results = results | records
        d += 1
    return results


def _df_edges_from_parents(
    record: Registry, field: str, distance: int, children: bool = False
):
    """Construct a DataFrame of edges as the input of graphviz.Digraph."""
    key = "children" if children else "parents"
    parents = _get_parents(
        record=record, field=field, distance=distance, children=children
    )
    all = record.__class__.objects
    records = parents | all.filter(id=record.id)
    df = records.distinct().df(include=[f"{key}__id"])
    df_edges = df[[f"{key}__id"]]
    df_edges = df_edges.explode(f"{key}__id")
    df_edges.index.name = "target"
    df_edges = df_edges.reset_index()
    df_edges.dropna(axis=0, inplace=True)
    df_edges.rename(columns={f"{key}__id": "source"}, inplace=True)
    df_edges = df_edges.drop_duplicates()

    # colons messes with the node formatting:
    # https://graphviz.readthedocs.io/en/stable/node_ports.html
    df_edges["source_record"] = df_edges["source"].apply(lambda x: all.get(id=x))
    df_edges["target_record"] = df_edges["target"].apply(lambda x: all.get(id=x))
    if record.__class__.__name__ == "Transform":
        df_edges["source_label"] = df_edges["source_record"].apply(
            _label_file_run_transform
        )
        df_edges["target_label"] = df_edges["target_record"].apply(
            _label_file_run_transform
        )
    else:
        df_edges["source_label"] = df_edges["source_record"].apply(
            lambda x: x.__getattribute__(field)
        )
        df_edges["target_label"] = df_edges["target_record"].apply(
            lambda x: x.__getattribute__(field)
        )
    return df_edges


def _add_emoji(record: Registry, label: str):
    if record.__class__.__name__ == "Transform":
        emoji = TRANSFORM_EMOJIS.get(record.type, "💫")
    elif record.__class__.__name__ == "Run":
        emoji = TRANSFORM_EMOJIS.get(record.transform.type, "💫")
    else:
        emoji = ""
    return f"{emoji} {label}"


def _get_all_parent_runs(file: File):
    """Get all input file runs recursively."""
    all_runs = {file.run}

    runs = [file.run]
    while any([r.input_files.exists() for r in runs if r is not None]):
        inputs = []
        for r in runs:
            inputs += r.input_files.all()
        runs = [f.run for f in inputs]
        all_runs.update(runs)
    return all_runs


def _get_all_child_runs(file: File):
    """Get all output file runs recursively."""
    all_runs: Set[Run] = set()

    runs = {f.run for f in file.run.output_files.all()}
    while runs.difference(all_runs):
        all_runs.update(runs)
        child_runs: Set[Run] = set()
        for r in runs:
            child_runs.update(
                Run.filter(input_files__id__in=r.output_files.list("id")).list()
            )
        runs = child_runs
    return all_runs


def _label_file_run_transform(record: Union[File, Run, Transform]):
    if isinstance(record, File):
        if record.description is None:
            name = record.key
        elif record.description.startswith("See dataset "):
            dataset_id = record.description.replace("See dataset ", "")
            dataset = Dataset.filter(id=dataset_id).one()
            name = dataset.name.replace("&", "&amp;")
        else:
            name = record.description.replace("&", "&amp;")

        return (
            rf'<{name}<BR/><FONT COLOR="GREY" POINT-SIZE="10"'
            rf' FACE="Monospace">id={record.id}<BR/>suffix={record.suffix}</FONT>>'
        )
    elif isinstance(record, Run):
        name = f'{record.transform.name.replace("&", "&amp;")}'
        return (
            rf'<{TRANSFORM_EMOJIS.get(str(record.transform.type), "💫")} {name}<BR/><FONT COLOR="GREY" POINT-SIZE="10"'  # noqa
            rf' FACE="Monospace">id={record.id}<BR/>type={record.transform.type},'
            rf" user={record.created_by.name}<BR/>run_at={format_field_value(record.run_at)}</FONT>>"  # noqa
        )
    elif isinstance(record, Transform):
        name = f'{record.name.replace("&", "&amp;")}'
        return (
            rf'<{TRANSFORM_EMOJIS.get(str(record.type), "💫")} {name}<BR/><FONT COLOR="GREY" POINT-SIZE="10"'  # noqa
            rf' FACE="Monospace">id={record.id}<BR/>type={record.type},'
            rf" user={record.created_by.name}<BR/>updated_at={format_field_value(record.updated_at)}</FONT>>"  # noqa
        )


def _df_edges_from_runs(all_runs: List[Run]):
    import pandas as pd

    df_values = []
    for run in all_runs:
        if run is not None:
            if run.input_files.exists():
                df_values.append((run.input_files.list(), run))
            if run.output_files.exists():
                df_values.append((run, run.output_files.list()))
    df = pd.DataFrame(df_values, columns=["source_record", "target_record"])
    df = df.explode("source_record")
    df = df.explode("target_record")
    df = df.drop_duplicates()
    df["source"] = [i.id for i in df["source_record"]]
    df["target"] = [i.id for i in df["target_record"]]
    df["source_label"] = df["source_record"].apply(_label_file_run_transform)
    df["target_label"] = df["target_record"].apply(_label_file_run_transform)
    return df


METHOD_NAMES = [
    "view_parents",
]

if _TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(HasParents, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, HasParents, globals())
