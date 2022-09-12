from pathlib import Path
from urllib.request import urlretrieve


def file_fcs() -> Path:
    """Return fcs file example."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/example.fcs", "example.fcs"
    )
    return Path(filepath)


def file_jpg_paradisi05() -> Path:
    """Return jpg file example."""
    filepath, _ = urlretrieve(
        "https://upload.wikimedia.org/wikipedia/commons/2/28/Laminopathic_nuclei.jpg",
        "paradisi05_laminopathic_nuclei.jpg",
    )
    return Path(filepath)


def dir_scrnaseq_cellranger() -> Path:
    """Directory with exemplary scrnaseq cellranger output."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/cellranger_run_001.zip",
    )
    from zipfile import ZipFile

    with ZipFile(filepath, "r") as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall(path=".")

    return Path("cellranger_run_001")


def file_mouse_sc_lymph_node() -> Path:
    """Mouse lymph node scRNA-seq dataset from EBI.

    Subsampled to 10k genes.

    From: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8414/
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/E-MTAB-8414.h5ad",
        "mouse_sc_lymph_node.h5ad",
    )
    return Path(filepath)


def schmidt22_crispra_gws_IFNG() -> Path:
    """CRISPRi screen dataset of Schmidt22.

    Originally from: https://zenodo.org/record/5784651
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/schmidt22-crispra-gws-IFNG.csv",
        "schmidt22-crispra-gws-IFNG.csv",
    )
    return Path(filepath)


def schmidt22_perturbseq() -> Path:
    """Perturb-seq dataset of Schmidt22.

    Subsampled and converted to h5ad from R file: https://zenodo.org/record/5784651

    To reproduce the subsample:
    >>> adata = sc.read('HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad')
    >>> adata.obs = adata.obs[['cluster_name']]
    >>> del adata.obsp
    >>> del adata.var['features']
    >>> del adata.obsm['X_pca']
    >>> del adata.uns
    >>> del adata.raw
    >>> del adata.varm
    >>> adata.obs = adata.obs.reset_index()
    >>> del adata.obs['index']
    >>> sc.pp.subsample(adata, 0.03)
    >>> adata.write('schmidt22_perturbseq.h5ad')
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/schmidt22_perturbseq.h5ad",
        "schmidt22_perturbseq.h5ad",
    )
    return Path(filepath)


def dir_scrnaseq_cellranger_schmidt22() -> Path:
    """BFXoutput directory containing Schmidt22_perturbseq."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/scrnaseq-cellranger-schmidt22.zip",
    )
    from zipfile import ZipFile

    with ZipFile(filepath, "r") as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall(path=".")

    return Path("scrnaseq-cellranger-schmidt22")
