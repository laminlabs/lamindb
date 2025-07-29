from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ... import Schema


def anndata_ensembl_gene_ids_and_valid_features_in_obs() -> Schema:
    """An `AnnData` schema validating Ensembl gene IDs and valid features in obs.

    .. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py
        :language: python
    """
    from ... import Schema

    docs_path = Path(__file__).parent.parent.parent.parent / "docs" / "scripts"
    if str(docs_path) not in sys.path:
        sys.path.append(str(docs_path))

    try:
        return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")
    except Schema.DoesNotExist:
        import define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs  # noqa

        try:
            return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")
        except Schema.DoesNotExist:
            importlib.reload(
                define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs
            )
            return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")
