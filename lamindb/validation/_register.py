from typing import Dict, List, Optional, Tuple, Union

import anndata as ad
import pandas as pd
from lamin_utils import colors, logger
from lnschema_core.types import FieldAttr

import lamindb as ln

from ._validate import (
    check_registry_organism,
    get_registry_instance,
    standardize_and_inspect,
)


def register_artifact(
    data: Union[pd.DataFrame, ad.AnnData],
    description: str,
    fields: Dict[str, FieldAttr],
    feature_field: FieldAttr,
    **kwargs,
) -> ln.Artifact:
    """Register all metadata with an Artifact.

    Args:
        data: The DataFrame or AnnData object to register.
        description: A description of the artifact.
        fields: A dictionary mapping obs_column to registry_field.
        feature_field: The registry field to validate variables index against.
        kwargs: Additional keyword arguments to pass to the registry model.

    Returns:
        The registered Artifact.
    """
    if isinstance(data, ad.AnnData):
        artifact = ln.Artifact.from_anndata(data, description=description)
        artifact.n_observations = data.n_obs
    elif isinstance(data, pd.DataFrame):
        artifact = ln.Artifact.from_df(data, description=description)
    else:
        raise ValueError("data must be a DataFrame or AnnData object")
    artifact.save()

    feature_kwargs: Dict = {}
    organism = check_registry_organism(
        feature_field.field.model, kwargs.pop("organism", None)
    )
    if organism is not None:
        feature_kwargs["organism"] = organism

    if isinstance(data, ad.AnnData):
        artifact.features.add_from_anndata(var_field=feature_field, **feature_kwargs)
    else:
        artifact.features.add_from_df(field=feature_field, **feature_kwargs)

    features = ln.Feature.lookup().dict()
    for feature_name, field in fields.items():
        feature = features.get(feature_name)
        registry = field.field.model
        filter_kwargs = kwargs.copy()
        organism = check_registry_organism(registry, organism)
        if organism is not None:
            filter_kwargs["organism"] = organism
        df = data.obs if isinstance(data, ad.AnnData) else data
        labels = registry.from_values(df[feature_name], field=field, **filter_kwargs)
        artifact.labels.add(labels, feature)

    slug = ln.setup.settings.instance.slug
    logger.success(f"registered artifact in {colors.italic(slug)}")
    if ln.setup.settings.instance.is_remote:
        logger.info(f"🔗 https://lamin.ai/{slug}/artifact/{artifact.uid}")

    return artifact


def register_labels(
    values: List[str],
    field: FieldAttr,
    feature_name: str,
    using: Optional[str] = None,
    validated_only: bool = True,
    kwargs: Optional[Dict] = None,
    df: Optional[pd.DataFrame] = None,
) -> None:
    """Register features or labels records in the default instance from the using instance.

    Args:
        values: A list of values to be registered as labels.
        field: The FieldAttr object representing the field for which labels are being registered.
        feature_name: The name of the feature to register.
        using: The name of the instance from which to transfer labels (if applicable).
        validated_only: If True, only register validated labels.
        kwargs: Additional keyword arguments to pass to the registry model.
        df: A DataFrame to register labels from.
    """
    filter_kwargs = {} if kwargs is None else kwargs.copy()
    registry = field.field.model
    if registry == ln.ULabel:
        validated_only = False

    organism = check_registry_organism(registry, filter_kwargs.pop("organism", None))
    if organism is not None:
        filter_kwargs["organism"] = organism

    verbosity = ln.settings.verbosity
    try:
        ln.settings.verbosity = "error"
        inspect_result_current = standardize_and_inspect(
            values=values, field=field, registry=registry, **filter_kwargs
        )
        if not inspect_result_current.non_validated:
            ln.settings.verbosity = verbosity
            return

        labels_registered: Dict = {"from public": [], "without reference": []}

        (
            labels_registered[f"from {using}"],
            non_validated_labels,
        ) = register_labels_from_using_instance(
            inspect_result_current.non_validated,
            field=field,
            using=using,
            kwargs=filter_kwargs,
        )

        public_records = (
            registry.from_values(non_validated_labels, field=field, **filter_kwargs)
            if non_validated_labels
            else []
        )
        ln.save(public_records)
        labels_registered["from public"] = [
            getattr(r, field.field.name) for r in public_records
        ]
        labels_registered["without reference"] = [
            i for i in non_validated_labels if i not in labels_registered["from public"]
        ]

        if not validated_only:
            non_validated_records = []
            if df is not None and registry == ln.Feature:
                non_validated_records = ln.Feature.from_df(df)
            else:
                if "organism" in filter_kwargs:
                    filter_kwargs["organism"] = _register_organism(name=organism)
                for value in labels_registered["without reference"]:
                    filter_kwargs[field.field.name] = value
                    if registry == ln.Feature:
                        filter_kwargs["type"] = "category"
                    non_validated_records.append(registry(**filter_kwargs))
            ln.save(non_validated_records)

        if registry == ln.ULabel and field.field.name == "name":
            register_ulabels_with_parent(values, field=field, feature_name=feature_name)
    finally:
        ln.settings.verbosity = verbosity

    log_registered_labels(
        labels_registered,
        feature_name=feature_name,
        model_field=f"{registry.__name__}.{field.field.name}",
        validated_only=validated_only,
    )


def log_registered_labels(
    labels_registered: Dict,
    feature_name: str,
    model_field: str,
    validated_only: bool = True,
) -> None:
    """Log the registered labels."""
    labels_type = "features" if feature_name == "feature" else "labels"
    model_field = colors.italic(model_field)
    for key, labels in labels_registered.items():
        if not labels:
            continue

        if key == "without reference" and validated_only:
            msg = colors.yellow(
                f"{len(labels)} non-validated {labels_type} are not registered with {model_field}: {labels}!"
            )
            lookup_print = f".lookup().['{feature_name}']"
            msg += f"\n      → to lookup categories, use {lookup_print}"
            msg += (
                f"\n      → to register, run {colors.yellow('register_features(validated_only=False)')}"
                if labels_type == "features"
                else f"\n      → to register, set {colors.yellow('validated_only=False')}"
            )
            logger.warning(msg)
        else:
            key = "" if key == "without reference" else f"{colors.green(key)} "
            logger.success(
                f"registered {len(labels)} {labels_type} {key}with {model_field}: {labels}"
            )


def register_ulabels_with_parent(
    values: List[str], field: FieldAttr, feature_name: str
) -> None:
    """Register a parent label for the given labels."""
    registry = field.field.model
    assert registry == ln.ULabel
    all_records = registry.from_values(values, field=field)
    is_feature = registry.filter(name=f"is_{feature_name}").one_or_none()
    if is_feature is None:
        is_feature = registry(name=f"is_{feature_name}")
        is_feature.save()
    is_feature.children.add(*all_records)


def register_labels_from_using_instance(
    values: List[str],
    field: FieldAttr,
    using: Optional[str] = None,
    kwargs: Optional[Dict] = None,
) -> Tuple[List[str], List[str]]:
    """Register features or labels records from the using instance.

    Args:
        values: A list of values to be registered as labels.
        field: The FieldAttr object representing the field for which labels are being registered.
        using: The name of the instance from which to transfer labels (if applicable).
        kwargs: Additional keyword arguments to pass to the registry model.

    Returns:
        A tuple containing the list of registered labels and the list of non-registered labels.
    """
    kwargs = kwargs or {}
    labels_registered = []
    not_registered = values

    if using is not None and using != "default":
        registry = field.field.model
        registry_using = get_registry_instance(registry, using)
        inspect_result_using = standardize_and_inspect(
            values=values, field=field, registry=registry_using, **kwargs
        )
        labels_using = registry_using.filter(
            **{f"{field.field.name}__in": inspect_result_using.validated}
        ).all()
        for label_using in labels_using:
            label_using.save()
            labels_registered.append(getattr(label_using, field.field.name))
        not_registered = inspect_result_using.non_validated

    return labels_registered, not_registered


def _register_organism(name: str):
    """Register an organism record."""
    import bionty as bt

    organism = bt.Organism.filter(name=name).one_or_none()
    if organism is None:
        organism = bt.Organism.from_public(name=name)
        if organism is None:
            raise ValueError(
                f"Organism '{name}' not found\n"
                f"      → please register it: bt.Organism(name='{name}').save()"
            )
        organism.save()
    return organism
