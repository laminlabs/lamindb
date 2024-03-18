from typing import Dict, List, Optional, Union

import anndata as ad
import pandas as pd
from lamin_utils import colors, logger
from lnschema_core.types import FieldAttr

import lamindb as ln

from ._validate import _registry_using, check_if_registry_needs_organism


def register_artifact(
    data: Union[pd.DataFrame, ad.AnnData],
    description: str,
    fields: Dict[str, FieldAttr],
    var_field: Optional[FieldAttr] = None,
    **kwargs,
):
    """Registers all metadata with an Artifact.

    Args:
        data: The DataFrame or AnnData object to register.
        description: A description of the artifact.
        fields: A dictionary mapping obs_column to registry_field.
        var_field: The registry field to validate variables index against.
        kwargs: Additional keyword arguments to pass to the registry model.
    """
    if isinstance(data, ad.AnnData):
        artifact = ln.Artifact.from_anndata(data, description=description)
        artifact.n_observations = data.n_obs
    elif isinstance(data, pd.DataFrame):
        artifact = ln.Artifact.from_df(data, description=description)
    else:
        raise ValueError("data must be a DataFrame or AnnData object")
    artifact.save()

    organism = kwargs.pop("organism", None)

    if isinstance(data, ad.AnnData):
        artifact.features.add_from_anndata(var_field=var_field, organism=organism)
    else:
        artifact.features.add_from_df()

    # link validated obs metadata
    features = ln.Feature.lookup().dict()
    for feature_name, field in fields.items():
        feature = features.get(feature_name)
        registry = field.field.model
        filter_kwargs = kwargs.copy()
        if check_if_registry_needs_organism(registry, organism):
            filter_kwargs["organism"] = organism
        df = data.obs if isinstance(data, ad.AnnData) else data
        labels = registry.from_values(df[feature_name], field=field, **filter_kwargs)
        artifact.labels.add(labels, feature)

    logger.print("\n\n🎉 registered artifact in LaminDB!\n")
    if ln.setup.settings.instance.is_remote:
        logger.print(
            f"🔗 https://lamin.ai/{ln.setup.settings.instance.slug}/artifact/{artifact.uid}"
        )

    return artifact


def register_labels(
    values: List[str],
    field: FieldAttr,
    feature_name: str,
    using: Optional[str] = None,
    validated_only: bool = True,
    kwargs: Dict = None,
):
    """Register features or labels records in the default instance from the using instance.

    Args:
        values: A list of values to be registered as labels.
        field: The FieldAttr object representing the field for which labels are being registered.
        feature_name: The name of the feature to register.
        using: The name of the instance from which to transfer labels (if applicable).
        validated_only: If True, only register validated labels.
        kwargs: Additional keyword arguments to pass to the registry model.
    """
    if kwargs is None:
        kwargs = {}
    registry = field.field.model

    check_if_registry_needs_organism(registry, kwargs.get("organism"))
    verbosity = ln.settings.verbosity
    try:
        ln.settings.verbosity = "error"
        # for labels that are registered in the using instance, transfer them to the current instance
        # first inspect the current instance
        inspect_result_current = registry.inspect(
            values, field=field, mute=True, **kwargs
        )
        if len(inspect_result_current.non_validated) == 0:
            # everything is validated in the current instance, no need to register
            ln.settings.verbosity = verbosity
            return

        labels_registered: Dict = {"from public": [], "without reference": []}

        # register labels from the using instance
        (
            labels_registered[f"from {using}"],
            non_validated_labels,
        ) = register_labels_from_using_instance(
            inspect_result_current.non_validated,
            field=field,
            using=using,
            kwargs=kwargs,
        )

        # for labels that are not registered in the using instance, register them in the current instance
        from_values_records = (
            registry.from_values(non_validated_labels, field=field, **kwargs)
            if len(non_validated_labels) > 0
            else []
        )
        ln.save(from_values_records)
        labels_registered["from public"] = [
            getattr(r, field.field.name) for r in from_values_records
        ]
        labels_registered["without reference"] = [
            i for i in non_validated_labels if i not in labels_registered["from public"]
        ]
        if not validated_only:
            non_validated_records = []
            for value in labels_registered["without reference"]:
                kwargs[field.field.name] = value
                if registry.__name__ == "Feature":
                    kwargs["type"] = "category"
                # register non-validated labels
                non_validated_records.append(registry(**kwargs))
            ln.save(non_validated_records)

        # for ulabels, also register a parent label: is_{feature_name}
        if registry == ln.ULabel and field.field.name == "name":
            register_ulabels_with_parent(values, field)
    finally:
        ln.settings.verbosity = verbosity
    log_registered_labels(
        labels_registered, feature_name=feature_name, validated_only=validated_only
    )


def log_registered_labels(
    labels_registered: Dict, feature_name: str, validated_only: bool = True
):
    """Log the registered labels."""
    for key, labels in labels_registered.items():
        if len(labels) > 0:
            if key == "without reference" and validated_only:
                msg = (
                    f"{len(labels)} non-validated labels are not registered: {labels}!\n"
                    "      → to lookup categories, use `.lookup().{feature_name}`\n"
                    "      → to register, set `validated_only=False`"
                )
                logger.warning(colors.yellow(msg))
                continue
            logger.success(
                f"registered {len(labels)} records {colors.green(key)}: {labels}"
            )


def register_ulabels_with_parent(values: List[str], field: FieldAttr):
    """Register a parent label for the given labels."""
    registry = field.field.model
    assert registry == ln.ULabel
    all_records = registry.from_values(values, field=field)
    is_feature = registry.filter(name=f"is_{field.field.name}").one_or_none()
    if is_feature is None:
        is_feature = registry(name=f"is_{field.field.name}")
        is_feature.save()
    # link all labels to the parent label
    is_feature.children.add(*all_records)


def register_labels_from_using_instance(
    values: List[str],
    field: FieldAttr,
    using: Optional[str] = None,
    kwargs: Dict = None,
):
    """Register features or labels records from the using instance.

    Args:
        values: A list of values to be registered as labels.
        field: The FieldAttr object representing the field for which labels are being registered.
        using: The name of the instance from which to transfer labels (if applicable).
        kwargs: Additional keyword arguments to pass to the registry model.
    """
    if kwargs is None:
        kwargs = {}
    labels_registered = []
    not_registered = values
    if using is not None and using != "default":
        registry = field.field.model
        registry_using = _registry_using(registry, using)
        # then inspect the using instance
        inspect_result_using = registry_using.inspect(
            values, field=field, mute=True, **kwargs
        )
        # register the labels that are validated in the using instance
        # TODO: filter kwargs
        labels_using = registry_using.filter(
            **{f"{field.field.name}__in": inspect_result_using.validated}
        ).all()
        for label_using in labels_using:
            label_using.save()
            labels_registered.append(getattr(label_using, field.field.name))
        not_registered = inspect_result_using.non_validated
    return labels_registered, not_registered
