from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

from django.db import connections
from lamin_utils import logger
from rich.table import Column, Table
from rich.text import Text
from rich.tree import Tree

from lamindb.models import CanCurate, Feature
from lamindb.models._from_values import _format_values
from lamindb.models.save import save
from lamindb.models.sqlrecord import (
    REGISTRY_UNIQUE_FIELD,
    get_name_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)

from ._describe import (
    NAME_WIDTH,
    TYPE_WIDTH,
    VALUES_WIDTH,
    describe_header,
    format_rich_tree,
)
from ._django import get_artifact_with_related, get_related_model
from ._relations import dict_related_model_to_related_name

if TYPE_CHECKING:
    from lamindb.models import Artifact, Collection, SQLRecord
    from lamindb.models.query_set import QuerySet

# we do not want to show records because this is a breaking change until all instances are migrated
# TODO: remove records from below once all instances are migrated
EXCLUDE_LABELS = {"feature_sets", "records"}


def _get_labels(
    obj, links: bool = False, instance: str | None = None
) -> dict[str, QuerySet]:
    """Get all labels associated with an object as a dictionary.

    This is a generic approach that uses django orm.
    """
    if obj.id is None:
        return {}

    labels = {}
    related_models = dict_related_model_to_related_name(
        obj.__class__, links=links, instance=instance
    )

    for _, related_name in related_models.items():
        if related_name not in EXCLUDE_LABELS and not related_name.startswith("_"):
            labels[related_name] = getattr(obj, related_name).all()
    return labels


def _get_labels_postgres(
    self: Artifact | Collection, m2m_data: dict | None = None
) -> dict[str, dict[int, str]]:
    """Get all labels associated with an artifact or collection as a dictionary.

    This is a postgres-specific approach that uses django Subquery.
    """
    if m2m_data is None:
        artifact_meta = get_artifact_with_related(self, include_m2m=True)
        m2m_data = artifact_meta.get("related_data", {}).get("m2m", {})
    return m2m_data


def describe_labels(
    self: Artifact | Collection,
    labels_data: dict | None = None,
    tree: Tree | None = None,
    as_subtree: bool = False,
):
    """Describe labels associated with an artifact or collection."""
    if not self._state.adding and connections[self._state.db].vendor == "postgresql":
        labels_data = _get_labels_postgres(self, labels_data)
    if not labels_data:
        labels_data = _get_labels(self, instance=self._state.db)

    # initialize tree
    if tree is None:
        tree = describe_header(self)
    if not labels_data:
        return tree

    labels_table = Table(
        Column("", style="", no_wrap=True, width=NAME_WIDTH),
        Column("", style="dim", no_wrap=True, width=TYPE_WIDTH),
        Column("", width=VALUES_WIDTH, no_wrap=True),
        show_header=False,
        box=None,
        pad_edge=False,
    )
    for related_name, labels in labels_data.items():
        if not labels or related_name == "feature_sets":
            continue
        if isinstance(labels, dict):  # postgres, labels are a dict[id, name]
            print_values = _format_values(labels.values(), n=10, quotes=False)
        else:  # labels are a QuerySet
            field = get_name_field(labels)
            print_values = _format_values(
                labels.values_list(field, flat=True), n=10, quotes=False
            )
        if print_values:
            related_model = get_related_model(self, related_name)
            type_str = related_model.__get_name_with_module__()
            labels_table.add_row(
                f".{related_name}", Text(type_str, style="dim"), print_values
            )

    labels_header = Text("Labels", style="bold green_yellow")
    if as_subtree:
        if labels_table.rows:
            labels_tree = Tree(labels_header, guide_style="dim")
            labels_tree.add(labels_table)
            return labels_tree
    else:
        if labels_table.rows:
            labels_tree = tree.add(labels_header)
            labels_tree.add(labels_table)
        return tree


def _save_validated_records(
    labels: QuerySet | list | dict,
) -> list[str]:
    if not labels:
        return []
    registry = labels[0].__class__
    field = (
        REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
        if not hasattr(registry, "_ontology_id_field")
        else registry._ontology_id_field
    )
    # if the field value is None, use uid field
    label_uids = [getattr(label, field) for label in labels if label is not None]
    # save labels from ontology_ids
    if hasattr(registry, "_ontology_id_field") and label_uids:
        try:
            records = registry.from_values(label_uids, field=field, mute=True)
            save([r for r in records if r._state.adding])
        except Exception:  # noqa: S110
            pass
        field = "uid"
        label_uids = [label.uid for label in labels if label is not None]

    if issubclass(registry, CanCurate):
        validated = registry.validate(label_uids, field=field, mute=True)
        new_labels = [
            label for label, is_valid in zip(labels, validated) if not is_valid
        ]
        return new_labels
    return list(labels)


def save_validated_records(
    labels: QuerySet | list | dict,
) -> list[str] | dict[str, list[str]]:
    """Save validated labels from public based on ontology_id_fields."""
    if isinstance(labels, dict):
        return {
            registry: _save_validated_records(registry_labels)
            for registry, registry_labels in labels.items()
        }
    return _save_validated_records(labels)


class LabelManager:
    """Label manager.

    This allows to manage untyped labels :class:`~lamindb.ULabel` and arbitrary
    typed labels (e.g., :class:`~bionty.CellLine`) and associate labels
    with features.
    """

    def __init__(self, host: Artifact | Collection) -> None:
        self._host = host

    def __repr__(self) -> str:
        return self.describe(return_str=True)

    def describe(self, return_str=True) -> str:
        """Describe the labels."""
        tree = describe_labels(self._host)
        return format_rich_tree(
            tree, fallback="no linked labels", return_str=return_str
        )

    def add(
        self,
        records: SQLRecord | list[SQLRecord] | QuerySet,
        feature: Feature | None = None,
    ) -> None:
        """Add one or several labels and associate them with a feature.

        Args:
            records: Label records to add.
            feature: Feature under which to group the labels.
        """
        from .artifact import add_labels

        return add_labels(self._host, records=records, feature=feature)

    def get(
        self,
        feature: Feature,
        mute: bool = False,
        flat_names: bool = False,
    ) -> QuerySet | dict[str, QuerySet] | list:
        """Get labels given a feature.

        Args:
            feature: Feature under which labels are grouped.
            mute: Show no logging.
            flat_names: Flatten list to names rather than returning records.
        """
        from .artifact import get_labels

        return get_labels(self._host, feature=feature, mute=mute, flat_names=flat_names)

    def add_from(self, data: Artifact | Collection, transfer_logs: dict = None) -> None:
        """Add labels from an artifact or collection to another artifact or collection.

        Examples:
            >>> artifact1 = ln.Artifact(pd.DataFrame(index=[0, 1])).save()
            >>> artifact2 = ln.Artifact(pd.DataFrame(index=[2, 3])).save()
            >>> ulabels = ln.ULabel.from_values(["Label1", "Label2"], field="name")
            >>> ln.save(ulabels)
            >>> labels = ln.ULabel.filter(name__icontains = "label").all()
            >>> artifact1.ulabels.set(labels)
            >>> artifact2.labels.add_from(artifact1)
        """
        if transfer_logs is None:
            transfer_logs = {"mapped": [], "transferred": [], "run": None}
        from lamindb import settings

        using_key = settings._using_key
        for related_name, labels in _get_labels(data, instance=data._state.db).items():
            labels = labels.all()
            if not labels.exists():
                continue
            # look for features
            data_name_lower = data.__class__.__name__.lower()
            labels_by_features: dict = defaultdict(list)
            features = set()
            new_labels = save_validated_records(labels)
            if len(new_labels) > 0:
                transfer_fk_to_default_db_bulk(
                    new_labels, using_key, transfer_logs=transfer_logs
                )
            for label in labels:
                keys: list = []
                # if the link table doesn't follow this convention, we'll ignore it
                if not hasattr(label, f"links_{data_name_lower}"):
                    key = None
                    keys.append(key)
                else:
                    links = (
                        getattr(label, f"links_{data_name_lower}")
                        .filter(**{f"{data_name_lower}_id": data.id})
                        .all()
                    )
                    for link in links:
                        if link.feature is not None:
                            features.add(link.feature)
                            key = link.feature.name
                        else:
                            key = None
                        keys.append(key)
                label_returned = transfer_to_default_db(
                    label,
                    using_key,
                    transfer_logs=transfer_logs,
                    transfer_fk=False,
                    save=True,
                )
                # TODO: refactor return value of transfer to default db
                if label_returned is not None:
                    label = label_returned
                for key in keys:
                    labels_by_features[key].append(label)
            # treat features
            new_features = save_validated_records(list(features))
            if len(new_features) > 0:
                transfer_fk_to_default_db_bulk(
                    new_features, using_key, transfer_logs=transfer_logs
                )
                for feature in new_features:
                    transfer_to_default_db(
                        feature,  # type: ignore
                        using_key,
                        transfer_logs=transfer_logs,
                        transfer_fk=False,
                    )
                save(new_features)  # type: ignore
            if hasattr(self._host, related_name):
                for feature_name, feature_labels in labels_by_features.items():
                    if feature_name is not None:
                        feature_id = Feature.get(name=feature_name).id
                    else:
                        feature_id = None
                    getattr(self._host, related_name).add(
                        *feature_labels, through_defaults={"feature_id": feature_id}
                    )

    def make_external(self, label: SQLRecord) -> None:
        """Make a label external, aka dissociate label from internal features.

        Args:
            label: Label record to make external.
        """
        d = dict_related_model_to_related_name(self._host)
        registry = label.__class__
        related_name = d.get(registry.__get_name_with_module__())
        link_model = getattr(self._host, related_name).through
        link_records = link_model.filter(
            artifact_id=self._host.id, **{f"{registry.__name__.lower()}_id": label.id}
        )
        features = link_records.values_list("feature__name", flat=True).distinct()
        s = "s" if len(features) > 1 else ""
        link_records.update(feature_id=None, feature_ref_is_name=None)
        logger.warning(
            f'{registry.__name__} "{getattr(label, label._name_field)}" is no longer associated with the following feature{s}:\n'
            f"{list(features)}"
        )
