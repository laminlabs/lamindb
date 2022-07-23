from lamindb_schema.core import dobject

from ..dev.db._insert import insert
from ..dev.file import h5ad_to_anndata


def anndata_to_df(adata, obs_or_var):
    """Get a df from AnnData."""
    if obs_or_var == "var":
        df = adata.var
    elif obs_or_var == "obs":
        df = adata.obs
    else:
        raise KeyError("Please specify obs or var!")
    return df


class annotate:
    """Feature annotation."""

    @classmethod
    def genes(
        cls,
        dobject: dobject,
        species: str,
        column=None,
        obs_or_var=None,
        geneset_name: str = None,
    ):
        """Annotate genes."""
        filekey = f"{dobject.id}-{dobject.v}{dobject.file_suffix}"

        if dobject.file_suffix == ".h5ad":
            adata = h5ad_to_anndata(filekey)
            df = anndata_to_df(adata, obs_or_var=obs_or_var)

            # create a geneset row
            genes = df.index.unique().values if column is None else df[column].unique()
            geneset_id = insert.genes(
                genes=genes, geneset_name=geneset_name, species=species
            )

        else:
            raise NotImplementedError

        return geneset_id

    @classmethod
    def proteinset(cls):
        NotImplementedError

    @classmethod
    def biosample(cls):
        raise NotImplementedError

    @classmethod
    def readout_type(cls):
        raise NotImplementedError
