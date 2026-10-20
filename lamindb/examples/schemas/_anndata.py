from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

from lamindb_setup.core import deprecated

if TYPE_CHECKING:
    from ... import Schema


@deprecated("Define a custom schema instead.")
def anndata_ensembl_gene_ids_and_valid_features_in_obs() -> Schema:
    """An `AnnData` schema validating Ensembl gene IDs and valid features in obs.

    .. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py
        :language: python
    """
    from ... import Schema

    try:
        return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")
    except Schema.DoesNotExist:
        from . import define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs  # noqa

        try:
            return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")
        except Schema.DoesNotExist:
            importlib.reload(
                define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs
            )
            return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")


def anndata_ensembl_gene_ids_and_obs_feature_type() -> Schema:
    """An `AnnData` schema validating Ensembl gene IDs and a feature type in obs.

    .. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_obs_feature_type.py
        :language: python
    """
    from ... import Schema

    name = "anndata_ensembl_gene_ids_and_obs_feature_type"
    try:
        return Schema.get(name=name)
    except Schema.DoesNotExist:
        from . import (  # noqa
            define_schema_anndata_ensembl_gene_ids_and_obs_feature_type,
        )

        try:
            return Schema.get(name=name)
        except Schema.DoesNotExist:
            importlib.reload(
                define_schema_anndata_ensembl_gene_ids_and_obs_feature_type
            )
            return Schema.get(name=name)
