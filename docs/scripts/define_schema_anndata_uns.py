import lamindb as ln

from define_schema_df_metadata import study_metadata_schema

anndata_uns_schema = ln.Schema(
    otype="AnnData",
    slots={
        "uns:study_metadata": study_metadata_schema,
    },
).save()
