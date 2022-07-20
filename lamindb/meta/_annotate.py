from ..admin.db._insert import insert
from ..dev.file import h5ad_to_anndata
from ..schema import dobject


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

    @staticmethod
    def gene(dobject: dobject, column=None, obs_or_var=None):
        """Annotate genes."""
        filekey = f"{dobject.id}-{dobject.v}{dobject.file_suffix}"

        if dobject.file_suffix == ".h5ad":
            adata = h5ad_to_anndata(filekey)
            df = anndata_to_df(adata, obs_or_var=obs_or_var)
            genes = df.index.unique().values if column is None else df[column].unique()
            gene_ids: dict = {}
            for i in genes:
                gene_ids[i] = insert.gene(i)
        else:
            raise NotImplementedError

        return gene_ids
