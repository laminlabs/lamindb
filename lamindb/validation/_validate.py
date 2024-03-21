from typing import Dict, Iterable, Optional

import pandas as pd
from anndata import AnnData
from lamin_utils import colors, logger
from lnschema_core import Registry
from lnschema_core.types import FieldAttr

from lamindb._from_values import _print_values


def get_registry_instance(registry: Registry, using: Optional[str] = None) -> Registry:
    """Get a registry instance using a specific instance."""
    if using is not None and using != "default":
        return registry.using(using)
    return registry


def standardize_and_inspect(
    values: Iterable[str], field: FieldAttr, registry: Registry, **kwargs
):
    """Standardize and inspect values using a registry."""
    if hasattr(registry, "standardize"):
        values = registry.standardize(values, field=field, mute=True, **kwargs)
    return registry.inspect(values, field=field, mute=True, **kwargs)


def check_registry_organism(
    registry: Registry, organism: Optional[str] = None
) -> Optional[str]:
    """Check if a registry needs an organism and return the organism name."""
    if hasattr(registry, "organism_id"):
        import bionty as bt

        if organism is None and bt.settings.organism is None:
            raise ValueError(
                f"{registry.__name__} registry requires an organism!\n"
                "      → please pass an organism name via organism="
            )
        return organism or bt.settings.organism.name
    return None


def validate_categories(
    values: Iterable[str],
    field: FieldAttr,
    feature_name: str,
    using: Optional[str] = None,
    **kwargs,
) -> bool:
    """Validate ontology terms in a pandas series using LaminDB registries."""
    model_field = f"{field.field.model.__name__}.{field.field.name}"
    logger.indent = ""
    logger.info(
        f"inspecting '{colors.bold(feature_name)}' by {colors.italic(model_field)}"
    )
    logger.indent = "   "

    registry = field.field.model
    filter_kwargs = {}
    organism = check_registry_organism(registry, kwargs.get("organism"))
    if organism is not None:
        filter_kwargs["organism"] = organism

    # Inspect the default instance
    inspect_result = standardize_and_inspect(
        values=values, field=field, registry=registry, **filter_kwargs
    )
    non_validated = inspect_result.non_validated

    if using is not None and using != "default" and non_validated:
        registry = get_registry_instance(registry, using)
        # Inspect the using instance
        inspect_result = standardize_and_inspect(
            values=non_validated, field=field, registry=registry, **filter_kwargs
        )
        non_validated = inspect_result.non_validated

    n_non_validated = len(non_validated)
    if n_non_validated == 0:
        logger.success(f"all {feature_name}s are validated")
        return True
    else:
        are = "are" if n_non_validated > 1 else "is"
        print_values = _print_values(non_validated)
        feature_name_print = f".register_labels('{feature_name}')"
        warning_message = (
            f"{colors.yellow(f'{n_non_validated} terms')} {are} not validated: "
            f"{colors.yellow(print_values)}\n      → register terms via "
            f"{colors.yellow(feature_name_print)}"
        )
        logger.warning(warning_message)
        logger.indent = ""
        return False


def validate_categories_in_df(
    df: pd.DataFrame,
    fields: Dict[str, FieldAttr],
    using: Optional[str] = None,
    **kwargs,
) -> bool:
    """Validate categories in DataFrame columns using LaminDB registries."""
    validated = True
    for feature_name, field in fields.items():
        validated &= validate_categories(
            df[feature_name],
            field=field,
            feature_name=feature_name,
            using=using,
            **kwargs,
        )
    return validated


def validate_anndata(
    adata: AnnData,
    var_field: FieldAttr,
    obs_fields: Dict[str, FieldAttr],
    using: Optional[str] = None,
    **kwargs,
) -> bool:
    """Inspect metadata in an AnnData object using LaminDB registries."""
    if using is not None and using != "default":
        logger.important(
            f"validating metadata using registries of instance {colors.italic(using)}"
        )

    validated_var = validate_categories(
        adata.var.index,
        field=var_field,
        feature_name="variables",
        using=using,
        **kwargs,
    )
    validated_obs = validate_categories_in_df(
        adata.obs, fields=obs_fields, using=using, **kwargs
    )
    return validated_var and validated_obs
