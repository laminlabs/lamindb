from typing import Dict, List, Union

from lamin_utils import colors
from lnschema_core.models import Data, Dataset, Feature, File, Registry

from .._feature_set import dict_related_model_to_related_name
from .._from_values import _print_values
from .._query_set import QuerySet
from .._registry import get_default_str_field


def get_labels_as_dict(self: Data):
    labels = {}
    for related_model, related_name in dict_related_model_to_related_name(
        self.__class__
    ).items():
        if related_name in {"feature_sets", "files", "input_of", "datasets"}:
            continue
        labels[related_name] = (related_model, self.__getattribute__(related_name))
    return labels


def print_labels(self: Data):
    labels_msg = ""
    for related_name, (related_model, labels) in get_labels_as_dict(self).items():
        if labels.exists():
            n = labels.count()
            field = get_default_str_field(labels)
            print_values = _print_values(labels.list(field), n=10)
            labels_msg += f"  🏷️ {related_name} ({n}, {colors.italic(related_model)}): {print_values}\n"  # noqa
    if len(labels_msg) > 0:
        return f"{colors.green('Labels')}:\n{labels_msg}"
    else:
        return ""


class LabelManager:
    """Label manager (:attr:`~lamindb.dev.Data.labels`).

    This allows to manage untyped labels :class:`~lamindb.ULabel` and arbitrary
    typed labels (e.g., :class:`~lnschema_bionty.CellLine`) and associate labels
    with features.

    See :class:`~lamindb.dev.Data` for more information.
    """

    def __init__(self, host: Union[File, Dataset]):
        self._host = host

    def __repr__(self) -> str:
        msg = print_labels(self._host)
        if len(msg) > 0:
            return msg
        else:
            return "no linked labels"

    def add(
        self,
        records: Union[Registry, List[Registry], QuerySet],
        feature: Feature,
    ) -> None:
        """Add one or several labels and associate them with a feature.

        Args:
            records: Label records to add.
            feature: Feature under which to group the labels.
            field: Field to parse iterable with from_values.
        """
        from ._data import add_labels

        return add_labels(self._host, records=records, feature=feature)

    def get(
        self,
        feature: Feature,
        mute: bool = False,
        flat_names: bool = False,
    ) -> Union[QuerySet, Dict[str, QuerySet], List]:
        """Get labels given a feature.

        Args:
            feature: Feature under which labels are grouped.
            mute: Show no logging.
            flat_names: Flatten list to names rather than returning records.
        """
        from ._data import get_labels

        return get_labels(self._host, feature=feature, mute=mute, flat_names=flat_names)

    def add_from(self, data: Data):
        """Transfer labels from a file or dataset.

        Examples:
            >>> file1 = ln.File(pd.DataFrame(index=[0, 1]))
            >>> file1.save()
            >>> file2 = ln.File(pd.DataFrame(index=[2, 3]))
            >>> file2.save()
            >>> ulabels = ln.ULabel.from_values(["Label1", "Label2"], field="name")
            >>> ln.save(ulabels)
            >>> labels = ln.ULabel.filter(name__icontains = "label").all()
            >>> file1.ulabels.set(labels)
            >>> file2.labels.add_from(file1)
        """
        features_lookup = Feature.lookup()
        for _, feature_set in data.features._feature_set_by_slot.items():
            if feature_set.registry == "core.Feature":
                df_slot = feature_set.features.df()
                for _, row in df_slot.iterrows():
                    if row["type"] == "category" and row["registries"] is not None:
                        labels = data.labels.get(
                            getattr(features_lookup, row["name"]), mute=True
                        )
                        self._host.labels.add(
                            labels, feature=getattr(features_lookup, row["name"])
                        )
        # for now, have this be duplicated, need to disentangle above
        for related_name, (_, labels) in get_labels_as_dict(data).items():
            getattr(self._host, related_name).add(*labels.all())
