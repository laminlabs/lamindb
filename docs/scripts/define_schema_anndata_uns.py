import lamindb as ln

from define_schema_df_metadata import study_metadata_schema

anndata_uns_schema = ln.Schema(
    name="anndata_study_metadata",
    otype="AnnData",
    slots={
        "uns": study_metadata_schema,
    },
).save()
