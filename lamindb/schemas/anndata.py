"""AnnData schemas.

.. autosummary::
   :toctree: .

   ensembl_gene_ids_and_valid_features_in_obs

"""

from .. import Schema


def ensembl_gene_ids_and_valid_features_in_obs() -> Schema:
    """Return a schema for an AnnData with Ensembl gene IDs and valid features in obs.

    .. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py
        :language: python
    """
    import sys
    from pathlib import Path

    docs_path = Path(__file__).parent.parent.parent / "docs" / "scripts"
    if str(docs_path) not in sys.path:
        sys.path.append(str(docs_path))

    import define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs  # noqa

    return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")
