from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from django.db import connections
from lamin_utils import colors, logger
from lnschema_core.models import CanValidate, Feature

from lamindb._from_values import _print_values
from lamindb._record import (
    REGISTRY_UNIQUE_FIELD,
    get_name_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)
from lamindb._save import save

from ._django import get_artifact_with_related, get_related_model
from ._settings import settings
from .schema import dict_related_model_to_related_name

if TYPE_CHECKING:
    from lnschema_core.models import Artifact, Collection, Record

    from lamindb._query_set import QuerySet

LABELS_EXCLUDE_SET = {"feature_sets"}


def get_labels_as_dict(
    self: Artifact | Collection, links: bool = False, instance: str | None = None
) -> dict:
    labels = {}  # type: ignore
    if self.id is None:
        return labels
    for related_model_name, related_name in dict_related_model_to_related_name(
        self.__class__, links=links, instance=instance
    ).items():
        if related_name not in LABELS_EXCLUDE_SET and not related_name.startswith("_"):
            labels[related_name] = (
                related_model_name,
                getattr(self, related_name).all(),
            )
    return labels


def _print_labels_postgres(
    self: Artifact | Collection, m2m_data: dict | None = None, print_types: bool = False
) -> str:
    labels_msg = ""
    if not m2m_data:
        artifact_meta = get_artifact_with_related(self, include_m2m=True)
        m2m_data = artifact_meta.get("related_data", {}).get("m2m", {})
    if m2m_data:
        for related_name, labels in m2m_data.items():
            if not labels or related_name == "feature_sets":
                continue
            related_model = get_related_model(self, related_name)
            print_values = _print_values(labels.values(), n=10)
            type_str = f": {related_model}" if print_types else ""
            labels_msg += f"    .{related_name}{type_str} = {print_values}\n"
    return labels_msg


def print_labels(
    self: Artifact | Collection,
    m2m_data: dict | None = None,
    print_types: bool = False,
):
    if not self._state.adding and connections[self._state.db].vendor == "postgresql":
        labels_msg = _print_labels_postgres(self, m2m_data, print_types)
    else:
        labels_msg = ""
        for related_name, (related_model, labels) in get_labels_as_dict(
            self, instance=self._state.db
        ).items():
            field = get_name_field(labels)
            labels_list = list(labels.values_list(field, flat=True))
            if len(labels_list) > 0:
                print_values = _print_values(labels_list, n=10)
                type_str = f": {related_model}" if print_types else ""
                labels_msg += f"    .{related_name}{type_str} = {print_values}\n"

    msg = ""
    if labels_msg:
        msg += f"  {colors.italic('Labels')}\n"
        msg += labels_msg
    return msg


# Alex: is this a label transfer function?
def validate_labels(labels: QuerySet | list | dict):
    def validate_labels_registry(
        labels: QuerySet | list | dict,
    ) -> tuple[list[str], list[str]]:
        if len(labels) == 0:
            return [], []
        registry = labels[0].__class__
        field = REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
        if hasattr(registry, "_ontology_id_field"):
            field = registry._ontology_id_field
        # if the field value is None, use uid field
        label_uids = np.array(
            [getattr(label, field) for label in labels if label is not None]
        )
        # save labels from ontology_ids
        if hasattr(registry, "_ontology_id_field") and len(label_uids) > 0:
            try:
                labels_records = registry.from_values(label_uids, field=field)
                save([r for r in labels_records if r._state.adding])
            except Exception:  # noqa S110
                pass
            field = "uid"
            label_uids = np.array(
                [getattr(label, field) for label in labels if label is not None]
            )
        if issubclass(registry, CanValidate):
            validated = registry.validate(label_uids, field=field, mute=True)
            validated_uids = label_uids[validated]
            validated_labels = registry.filter(
                **{f"{field}__in": validated_uids}
            ).list()
            new_labels = [labels[int(i)] for i in np.argwhere(~validated).flatten()]
        else:
            validated_labels = []
            new_labels = list(labels)
        return validated_labels, new_labels

    if isinstance(labels, dict):
        result = {}
        for registry, labels_registry in labels.items():
            result[registry] = validate_labels_registry(labels_registry)
    else:
        return validate_labels_registry(labels)


class LabelManager:
    """Label manager.

    This allows to manage untyped labels :class:`~lamindb.ULabel` and arbitrary
    typed labels (e.g., :class:`~bionty.CellLine`) and associate labels
    with features.
    """

    def __init__(self, host: Artifact | Collection):
        self._host = host

    def __repr__(self) -> str:
        msg = print_labels(self._host)
        if len(msg) > 0:
            return msg
        else:
            return "no linked labels"

    def add(
        self,
        records: Record | list[Record] | QuerySet,
        feature: Feature | None = None,
    ) -> None:
        """Add one or several labels and associate them with a feature.

        Args:
            records: Label records to add.
            feature: Feature under which to group the labels.
        """
        from ._data import add_labels

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
        from ._data import get_labels

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
        using_key = settings._using_key
        for related_name, (_, labels) in get_labels_as_dict(
            data, instance=data._state.db
        ).items():
            labels = labels.all()
            if not labels.exists():
                continue
            # look for features
            data_name_lower = data.__class__.__name__.lower()
            labels_by_features = defaultdict(list)
            features = set()
            _, new_labels = validate_labels(labels)
            if len(new_labels) > 0:
                transfer_fk_to_default_db_bulk(
                    new_labels, using_key, transfer_logs=transfer_logs
                )
            for label in labels:
                # if the link table doesn't follow this convention, we'll ignore it
                if not hasattr(label, f"links_{data_name_lower}"):
                    key = None
                else:
                    link = getattr(label, f"links_{data_name_lower}").get(
                        **{f"{data_name_lower}_id": data.id}
                    )
                    if link.feature is not None:
                        features.add(link.feature)
                        key = link.feature.name
                    else:
                        key = None
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
                labels_by_features[key].append(label)
            # treat features
            _, new_features = validate_labels(list(features))
            if len(new_features) > 0:
                transfer_fk_to_default_db_bulk(
                    new_features, using_key, transfer_logs=transfer_logs
                )
                for feature in new_features:
                    transfer_to_default_db(
                        feature,
                        using_key,
                        transfer_logs=transfer_logs,
                        transfer_fk=False,
                    )
                save(new_features)
            if hasattr(self._host, related_name):
                for feature_name, labels in labels_by_features.items():
                    if feature_name is not None:
                        feature_id = Feature.get(name=feature_name).id
                    else:
                        feature_id = None
                    getattr(self._host, related_name).add(
                        *labels, through_defaults={"feature_id": feature_id}
                    )

    def make_external(self, label: Record):
        """Make a label external, aka dissociate label from internal features.

        Args:
            label: Label record to make external.
        """
        d = dict_related_model_to_related_name(self._host)
        registry = label.__class__
        related_name = d.get(registry.__get_name_with_schema__())
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
