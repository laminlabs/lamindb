from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Dict

import numpy as np
from lamin_utils import colors, logger
from lnschema_core.models import Feature

from lamindb._from_values import _print_values
from lamindb._registry import (
    REGISTRY_UNIQUE_FIELD,
    get_default_str_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)
from lamindb._save import save

from ._settings import settings
from .schema import dict_related_model_to_related_name

if TYPE_CHECKING:
    from lnschema_core.models import Artifact, Collection, Data, Registry

    from lamindb._query_set import QuerySet


def get_labels_as_dict(self: Data, links: bool = False):
    exclude_set = {
        "feature_sets",
        "unordered_artifacts",
        "input_of",
        "collections",
        "source_code_of",
        "report_of",
        "environment_of",
        "collection_links",
        "artifact_links",
        "feature_set_links",
        "previous_runs",
        "feature_values",
    }
    labels = {}  # type: ignore
    if self.id is None:
        return labels
    for related_model_name, related_name in dict_related_model_to_related_name(
        self.__class__, links=links
    ).items():
        if related_name not in exclude_set:
            labels[related_name] = (
                related_model_name,
                getattr(self, related_name).all(),
            )
    return labels


def print_labels(self: Data, field: str = "name", print_types: bool = False):
    labels_msg = ""
    for related_name, (related_model, labels) in get_labels_as_dict(self).items():
        # there is a try except block here to deal with schema inconsistencies
        # during transfer between instances
        try:
            labels_list = list(labels.values_list(field, flat=True))
            if len(labels_list) > 0:
                get_default_str_field(labels)
                print_values = _print_values(labels_list, n=10)
                type_str = f": {related_model}" if print_types else ""
                labels_msg += f"    .{related_name}{type_str} = {print_values}\n"
        except Exception:
            continue
    msg = ""
    if labels_msg:
        msg += f"  {colors.italic('Labels')}\n"
        msg += labels_msg
    return msg


def transfer_add_labels(
    labels, features_lookup_self, self, feature_name, parents: bool = True
):
    def transfer_single_registry(validated_labels, new_labels):
        # here the new labels are transferred to the self db
        if len(new_labels) > 0:
            transfer_fk_to_default_db_bulk(new_labels, using_key=None)
            for label in new_labels:
                transfer_to_default_db(
                    label, using_key=None, mute=True, transfer_fk=False
                )
            # not saving parents for Organism during transfer
            registry = new_labels[0].__class__
            logger.info(f"saving {len(new_labels)} new {registry.__name__} records")
            save(new_labels)
        # link labels records from self db
        self._host.labels.add(
            validated_labels + new_labels,
            feature=features_lookup_self.get(feature_name),
        )

    # validate labels on the default db
    result = validate_labels(labels, parents=parents)
    if isinstance(result, Dict):
        for _, (validated_labels, new_labels) in result.items():
            transfer_single_registry(validated_labels, new_labels)
    else:
        transfer_single_registry(*result)


# Alex: is this a label transfer function?
def validate_labels(labels: QuerySet | list | dict, parents: bool = True):
    def validate_labels_registry(
        labels: QuerySet | list | dict, parents: bool = True
    ) -> tuple[list[str], list[str]]:
        if len(labels) == 0:
            return [], []
        registry = labels[0].__class__
        field = REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
        if hasattr(registry, "ontology_id") and parents:
            field = "ontology_id"
        elif hasattr(registry, "ensembl_gene_id"):
            field = "ensembl_gene_id"
        elif hasattr(registry, "uniprotkb_id"):
            field = "uniprotkb_id"
        if registry.__get_name_with_schema__() == "bionty.Organism":
            parents = False
        # if the field value is None, use uid field
        label_uids = np.array(
            [getattr(label, field) for label in labels if label is not None]
        )
        # save labels from ontology_ids so that parents are populated
        if field == "ontology_id" and len(label_uids) > 0:
            try:
                records = registry.from_values(label_uids, field=field)
                if len(records) > 0:
                    save(records, parents=parents)
            except Exception:
                pass
            field = "uid"
            label_uids = np.array(
                [getattr(label, field) for label in labels if label is not None]
            )
        validated = registry.validate(label_uids, field=field, mute=True)
        validated_uids = label_uids[validated]
        validated_labels = registry.filter(**{f"{field}__in": validated_uids}).list()
        new_labels = [labels[int(i)] for i in np.argwhere(~validated).flatten()]
        return validated_labels, new_labels

    if isinstance(labels, Dict):
        result = {}
        for registry, labels_registry in labels.items():
            result[registry] = validate_labels_registry(
                labels_registry, parents=parents
            )
    else:
        return validate_labels_registry(labels, parents=parents)


class LabelManager:
    """Label manager (:attr:`~lamindb.core.Data.labels`).

    This allows to manage untyped labels :class:`~lamindb.ULabel` and arbitrary
    typed labels (e.g., :class:`~bionty.CellLine`) and associate labels
    with features.

    See :class:`~lamindb.core.Data` for more information.
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
        records: Registry | list[Registry] | QuerySet,
        feature: Feature | None = None,
    ) -> None:
        """Add one or several labels and associate them with a feature.

        Args:
            records: Label records to add.
            feature: Feature under which to group the labels.
        """
        from ._data import add_labels

        return add_labels(self._host, records=records, feature=feature)

    def remove(self, record: Registry):
        """Remove a label."""
        d = dict_related_model_to_related_name(self._host.__class__)
        related_name = d[record.__class__]
        getattr(self._host, related_name).remove(record)

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

    def add_from(self, data: Data, parents: bool = True) -> None:
        """Add labels from an artifact or collection to another artifact or collection.

        Examples:
            >>> file1 = ln.Artifact(pd.DataFrame(index=[0, 1]))
            >>> file1.save()
            >>> file2 = ln.Artifact(pd.DataFrame(index=[2, 3]))
            >>> file2.save()
            >>> ulabels = ln.ULabel.from_values(["Label1", "Label2"], field="name")
            >>> ln.save(ulabels)
            >>> labels = ln.ULabel.filter(name__icontains = "label").all()
            >>> file1.ulabels.set(labels)
            >>> file2.labels.add_from(file1)
        """
        from django.db.utils import ProgrammingError

        using_key = settings._using_key
        for related_name, (_, labels) in get_labels_as_dict(data).items():
            labels = labels.all()
            try:
                if not labels.exists():
                    continue
                # look for features
                data_name_lower = data.__class__.__name__.lower()
                labels_by_features = defaultdict(list)
                features = set()
                _, new_labels = validate_labels(labels, parents=parents)
                if len(new_labels) > 0:
                    transfer_fk_to_default_db_bulk(new_labels, using_key)
                for label in labels:
                    # if the link table doesn't follow this convention, we'll ignore it
                    if not hasattr(label, f"{data_name_lower}_links"):
                        key = None
                    else:
                        link = getattr(label, f"{data_name_lower}_links").get(
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
                        mute=True,
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
                    transfer_fk_to_default_db_bulk(new_features, using_key)
                    for feature in new_features:
                        transfer_to_default_db(
                            feature, using_key, mute=True, transfer_fk=False
                        )
                    save(new_features, parents=parents)
                if hasattr(self._host, related_name):
                    for feature_name, labels in labels_by_features.items():
                        if feature_name is not None:
                            feature_id = Feature.filter(name=feature_name).one().id
                        else:
                            feature_id = None
                        getattr(self._host, related_name).add(
                            *labels, through_defaults={"feature_id": feature_id}
                        )
            # ProgrammingError is raised when schemas don't match between source and target instances
            except ProgrammingError:
                continue
