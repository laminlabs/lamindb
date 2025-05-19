from ... import Schema


def anndata_ensembl_gene_ids_and_valid_features_in_obs() -> Schema:
    """Return a schema for an AnnData with Ensembl gene IDs and valid features in obs.

    .. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py
        :language: python
    """
    import subprocess
    from pathlib import Path

    docs_path = Path(__file__).parent.parent.parent.parent / "docs" / "scripts"
    subprocess.run(
        [
            "python",
            str(
                docs_path
                / "define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py"
            ),
        ],
        check=True,
    )

    return Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs")
