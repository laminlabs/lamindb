from pathlib import Path
from typing import Union
from urllib.request import urlretrieve

import anndata as ad
import pandas as pd
from lnschema_core.dev import id


def file_fcs() -> Path:
    """Return fcs file example."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/example.fcs", "example.fcs"
    )
    return Path(filepath)


def file_jpg_paradisi05() -> Path:
    """Return jpg file example.

    Originally from: https://upload.wikimedia.org/wikipedia/commons/2/28/Laminopathic_nuclei.jpg
    """  # noqa
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/Laminopathic_nuclei.jpg",
        "paradisi05_laminopathic_nuclei.jpg",
    )
    return Path(filepath)


def file_fastq() -> Path:
    """Mini mock fastq file."""
    with open("./input.fastq.gz", "w") as f:
        f.write("Mock fastq file.")
    return Path("./input.fastq.gz")


def file_bam() -> Path:
    """Mini mock bam file."""
    with open("./output.bam", "w") as f:
        f.write("Mock bam file.")
    return Path("./output.bam")


def file_mini_csv() -> Path:
    """Mini csv file."""
    filename = Path("./mini.csv")
    df = pd.DataFrame([1, 2, 3], columns=["test"])
    df.to_csv(filename, index=False)
    return filename


def dir_scrnaseq_cellranger() -> Path:
    """Directory with exemplary scrnaseq cellranger input and output."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/cellranger_run_001.zip"
    )
    from zipfile import ZipFile

    with ZipFile(filepath, "r") as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall(path=".")

    return Path("cellranger_run_001")


def anndata_mouse_sc_lymph_node() -> ad.AnnData:
    """Mouse lymph node scRNA-seq dataset from EBI.

    Subsampled to 10k genes.

    From: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8414/
    """
    filepath, _ = urlretrieve("https://lamindb-test.s3.amazonaws.com/E-MTAB-8414.h5ad")
    return ad.read(filepath)


def anndata_pbmc68k_reduced() -> ad.AnnData:
    """Modified from scanpy.datasets.pbmc68k_reduced()."""
    import scanpy as sc

    pbmc68k = sc.datasets.pbmc68k_reduced()
    pbmc68k.obs.rename(columns={"bulk_labels": "cell_type"}, inplace=True)
    pbmc68k.obs["cell_type"] = pbmc68k.obs["cell_type"].cat.rename_categories(
        {"Dendritic": "Dendritic cells", "CD14+ Monocyte": "CD14+ Monocytes"}
    )
    del pbmc68k.obs["G2M_score"]
    del pbmc68k.obs["S_score"]
    del pbmc68k.obs["phase"]
    del pbmc68k.obs["n_counts"]
    del pbmc68k.var["dispersions"]
    del pbmc68k.var["dispersions_norm"]
    del pbmc68k.var["means"]
    del pbmc68k.uns["rank_genes_groups"]
    del pbmc68k.uns["bulk_labels_colors"]
    del pbmc68k.raw
    sc.pp.subsample(pbmc68k, fraction=0.1, random_state=123)
    return pbmc68k


def anndata_pbmc3k_processed() -> ad.AnnData:
    """Modified from scanpy.pbmc3k_processed()."""
    import scanpy as sc

    pbmc3k = sc.datasets.pbmc3k_processed()
    pbmc3k.obs.rename(columns={"louvain": "cell_type"}, inplace=True)
    return pbmc3k


def anndata_human_immune_cells() -> ad.AnnData:
    """Cross-tissue immune cell analysis reveals tissue-specific features in humans.

    From: https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3  # noqa
    Dataset: Global

    To reproduce the subsample:
    >>> adata = sc.read('Global.h5ad')
    >>> adata.obs = adata.obs[['donor_id', 'tissue', 'cell_type', 'assay', 'tissue_ontology_term_id', 'cell_type_ontology_term_id', 'assay_ontology_term_id']].copy() . # noqa
    >>> sc.pp.subsample(adata, fraction=0.005)
    >>> del adata.uns["development_stage_ontology_term_id_colors"]
    >>> del adata.uns["sex_ontology_term_id_colors"]
    >>> sc.write('human_immune.h5ad', adata)
    """
    filepath, _ = urlretrieve("https://lamindb-test.s3.amazonaws.com/human_immune.h5ad")
    return ad.read(filepath)


def generate_cell_ranger_files(
    sample_name: str, basedir: Union[str, Path] = "./", output_only: bool = True
):
    """Generate mock cell ranger outputs.

    Args:
        sample_name: name of the sample
        basedir: run directory
        output_only: only generate output files

    """
    basedir = Path(basedir)

    if not output_only:
        fastqdir = basedir / "fastq"
        fastqdir.mkdir(parents=True, exist_ok=True)
        fastqfile1 = fastqdir / f"{sample_name}_R1_001.fastq.gz"
        with open(fastqfile1, "w") as f:
            f.write(f"{id.base62(n_char=6)}")
        fastqfile2 = fastqdir / f"{sample_name}_R2_001.fastq.gz"
        fastqfile2.touch(exist_ok=True)
        with open(fastqfile2, "w") as f:
            f.write(f"{id.base62(n_char=6)}")

    sampledir = basedir / f"{sample_name}"
    for folder in ["raw_feature_bc_matrix", "filtered_feature_bc_matrix", "analysis"]:
        filedir = sampledir / folder
        filedir.mkdir(parents=True, exist_ok=True)

    for filename in [
        "web_summary.html",
        "metrics_summary.csv",
        "possorted_genome_bam.bam",
        "possorted_genome_bam.bam.bai",
        "molecule_info.h5",
        "cloupe.cloupe",
        "raw_feature_bc_matrix.h5",
        "raw_feature_bc_matrix/barcodes.tsv.gz",
        "raw_feature_bc_matrix/features.tsv.gz",
        "raw_feature_bc_matrix/matrix.mtx.gz",
        "filtered_feature_bc_matrix.h5",
        "filtered_feature_bc_matrix/barcodes.tsv.gz",
        "filtered_feature_bc_matrix/features.tsv.gz",
        "filtered_feature_bc_matrix/matrix.mtx.gz",
        "analysis/analysis.csv",
    ]:
        file = sampledir / filename
        with open(file, "w") as f:
            f.write(f"{id.base62(n_char=6)}")

    return sampledir


# def schmidt22_crispra_gws_IFNG() -> Path:
#     """CRISPRi screen dataset of Schmidt22.

#     Originally from: https://zenodo.org/record/5784651
#     """
#     filepath, _ = urlretrieve(
#         "https://lamindb-test.s3.amazonaws.com/schmidt22-crispra-gws-IFNG.csv",
#         "schmidt22-crispra-gws-IFNG.csv",
#     )
#     return Path(filepath)


# def schmidt22_perturbseq() -> Path:
#     """Perturb-seq dataset of Schmidt22.

#     Subsampled and converted to h5ad from R file: https://zenodo.org/record/5784651

#     To reproduce the subsample:
#     >>> adata = sc.read('HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad')
#     >>> adata.obs = adata.obs[['cluster_name']]
#     >>> del adata.obsp
#     >>> del adata.var['features']
#     >>> del adata.obsm['X_pca']
#     >>> del adata.uns
#     >>> del adata.raw
#     >>> del adata.varm
#     >>> adata.obs = adata.obs.reset_index()
#     >>> del adata.obs['index']
#     >>> sc.pp.subsample(adata, 0.03)
#     >>> adata.write('schmidt22_perturbseq.h5ad')
#     """
#     filepath, _ = urlretrieve(
#         "https://lamindb-test.s3.amazonaws.com/schmidt22_perturbseq.h5ad",
#         "schmidt22_perturbseq.h5ad",
#     )
#     return Path(filepath)


# def dir_scrnaseq_cellranger_schmidt22() -> Path:
#     """BFXoutput directory containing Schmidt22_perturbseq."""
#     filepath, _ = urlretrieve(
#         "https://lamindb-test.s3.amazonaws.com/scrnaseq-cellranger-schmidt22.zip",
#     )
#     from zipfile import ZipFile

#     with ZipFile(filepath, "r") as zipObj:
#         # Extract all the contents of zip file in current directory
#         zipObj.extractall(path=".")

#     return Path("scrnaseq-cellranger-schmidt22")
