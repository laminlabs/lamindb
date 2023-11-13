from pathlib import Path
from typing import Union
from urllib.request import urlretrieve

import anndata as ad
import numpy as np
import pandas as pd
from lnschema_core import ids
from upath import UPath

from .._settings import settings


def file_fcs() -> Path:
    """Example FCS file."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/example.fcs", "example.fcs"
    )
    return Path(filepath)


def file_fcs_alpert19(populate_registries: bool = False) -> Path:  # pragma: no cover
    """FCS file from Alpert19.

    Args:
        populate_registries: pre-populate metadata records to simulate existing registries  # noqa
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/Alpert19-070314-Mike-Study+15-2013-plate+1-15-004-1-13_cells_found.fcs",  # noqa
        "Alpert19.fcs",
    )
    if populate_registries:
        import lnschema_bionty as lb
        import readfcs

        import lamindb as ln

        verbosity = ln.settings.verbosity
        ln.settings.verbosity = "error"
        adata = readfcs.read(filepath)
        std = lb.CellMarker.bionty().standardize(adata.var.index)
        ln.save(
            lb.CellMarker.from_values(
                lb.CellMarker.bionty().inspect(std, "name").validated, "name"
            )
        )
        ln.Feature(
            name="assay", type="category", registries=[lb.ExperimentalFactor]
        ).save()
        ln.Feature(name="organism", type="category", registries=[lb.Organism]).save()
        ln.settings.verbosity = verbosity
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


def file_tsv_rnaseq_nfcore_salmon_merged_gene_counts(
    populate_registries: bool = False,
) -> Path:  # pragma: no cover
    """Gene counts table from nf-core RNA-seq pipeline.

    Output of: https://nf-co.re/rnaseq
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/salmon.merged.gene_counts.tsv",
        "salmon.merged.gene_counts.tsv",
    )
    if populate_registries:
        import lnschema_bionty as lb

        import lamindb as ln

        verbosity = ln.settings.verbosity
        ln.settings.verbosity = "error"
        ln.Feature(
            name="assay", type="category", registries=[lb.ExperimentalFactor]
        ).save()
        ln.Feature(name="organism", type="category", registries=[lb.Organism]).save()
        lb.ExperimentalFactor.from_bionty(ontology_id="EFO:0008896").save()
        ln.settings.verbosity = verbosity

    return Path(filepath)


def file_fastq(in_storage_root=False) -> Path:
    """Mini mock fastq file."""
    basedir = Path(".") if not in_storage_root else settings.storage
    filepath = basedir / "input.fastq.gz"
    with open(filepath, "w") as f:
        f.write("Mock fastq file.")
    return filepath


def file_bam(in_storage_root=False) -> Path:  # pragma: no cover
    """Mini mock bam file."""
    basedir = Path(".") if not in_storage_root else settings.storage
    filepath = basedir / "output.bam"
    with open(filepath, "w") as f:
        f.write("Mock bam file.")
    return filepath


def file_mini_csv(in_storage_root=False) -> Path:
    """Mini csv file."""
    basedir = Path(".") if not in_storage_root else settings.storage
    filepath = basedir / "mini.csv"
    df = pd.DataFrame([1, 2, 3], columns=["test"])
    df.to_csv(filepath, index=False)
    return filepath


def file_tiff_suo22():  # pragma: no cover
    """Image file from Suo22.

    Pair with anndata_suo22_Visium10X
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/F121_LP1_4LIV.tiff",
        "F121_LP1_4LIV.tiff",
    )
    Path("suo22/").mkdir(exist_ok=True)
    filepath = Path(filepath).rename("suo22/F121_LP1_4LIV.tiff")
    return Path(filepath)


def dir_iris_images() -> UPath:  # pragma: no cover
    """Directory with 3 studies of the Iris flower: 405 images & metadata.

    Based on: https://github.com/laminlabs/lamindb-dev-datasets/pull/2
    """
    return UPath("s3://lamindb-dev-datasets/iris_studies")


def anndata_mouse_sc_lymph_node(
    populate_registries: bool = False,
) -> ad.AnnData:  # pragma: no cover
    """Mouse lymph node scRNA-seq dataset from EBI.

    Subsampled to 10k genes.

    From: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8414/

    Args:
        populate_registries: pre-populate metadata records to simulate existing registries  # noqa
    """
    filepath, _ = urlretrieve("https://lamindb-test.s3.amazonaws.com/E-MTAB-8414.h5ad")
    adata = ad.read_h5ad(filepath)

    # The column names are a bit lengthy, let's abbreviate them:
    adata.obs.columns = (
        adata.obs.columns.str.replace("Sample Characteristic", "")
        .str.replace("Factor Value ", "Factor Value:", regex=True)
        .str.replace("Factor Value\[", "Factor Value:", regex=True)  # noqa
        .str.replace(" Ontology Term\[", "ontology_id:", regex=True)  # noqa
        .str.strip("[]")
        .str.replace("organism part", "tissue")
        .str.replace("organism", "organism")
        .str.replace("developmental stage", "developmental_stage")
        .str.replace("cell type", "cell_type")
        # the last one could be interesting, too
        # .str.replace("Factor Value:Ontology Term[inferred cell_type - authors labels", "cell_type_authors")  # noqa
    )
    # subset columns to only the ones with names
    columns = [
        col
        for col in adata.obs.columns
        if not col.startswith("ontology_id")
        and not col.startswith("Factor Value")
        and col != "strain"
    ]
    adata.obs = adata.obs[columns]

    # pre-populate registries
    if populate_registries:
        import lnschema_bionty as lb

        import lamindb as ln

        verbosity = ln.settings.verbosity
        ln.settings.verbosity = "error"
        auto_save_parents = lb.settings.auto_save_parents
        lb.settings.auto_save_parents = False
        # strain
        lb.ExperimentalFactor.from_bionty(ontology_id="EFO:0004472").save()
        # developmental stage
        lb.ExperimentalFactor.from_bionty(ontology_id="EFO:0001272").save()
        # tissue
        lb.Tissue.from_bionty(ontology_id="UBERON:0001542").save()
        # cell types
        ln.save(lb.CellType.from_values(["CL:0000115", "CL:0000738"], "ontology_id"))
        # assays
        ln.Feature(
            name="assay", type="category", registries=[lb.ExperimentalFactor]
        ).save()
        lb.ExperimentalFactor.from_bionty(ontology_id="EFO:0008913").save()
        # genes
        validated = lb.Gene.bionty(organism="mouse").validate(
            adata.var.index, field="ensembl_gene_id"
        )
        ln.save(
            lb.Gene.from_values(
                adata.var.index[validated][:-19],
                field="ensembl_gene_id",
                organism="mouse",
            )
        )
        # labels
        labels = []
        for col in ["sex", "age", "genotype", "immunophenotype"]:
            labels += [ln.ULabel(name=name) for name in adata.obs[col]]
        ln.save(labels)
        ln.settings.verbosity = verbosity
        lb.settings.auto_save_parents = auto_save_parents

    return adata


def anndata_pbmc68k_reduced() -> ad.AnnData:
    """Modified from scanpy.datasets.pbmc68k_reduced().

    This code was run::

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
        sc.pp.subsample(pbmc68k, fraction=0.1, random_state=123)
        pbmc68k.write("scrnaseq_pbmc68k_tiny.h5ad")
    """
    filepath, _ = urlretrieve(
        "https://lamindb-dev-datasets.s3.amazonaws.com/scrnaseq_pbmc68k_tiny.h5ad"
    )
    return ad.read_h5ad(filepath)


def anndata_file_pbmc68k_test() -> Path:
    """Modified from scanpy.datasets.pbmc68k_reduced().

    Additional slots were added for testing purposes. Returns the filepath.

    To reproduce::

        pbmc68k = ln.dev.datasets.anndata_pbmc68k_reduced()
        pbmc68k_test = pbmc68k[:30, :200].copy()
        pbmc68k_test.raw = pbmc68k_test[:, :100]
        pbmc68k_test.obsp["test"] = sparse.eye(pbmc68k_test.shape[0], format="csr")
        pbmc68k_test.varp["test"] = sparse.eye(pbmc68k_test.shape[1], format="csr")
        pbmc68k_test.layers["test"] = sparse.csr_matrix(pbmc68k_test.shape)
        pbmc68k_test.layers["test"][0] = 1.
        pbmc68k_test.write("pbmc68k_test.h5ad")
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/pbmc68k_test.h5ad", "pbmc68k_test.h5ad"
    )
    return Path(filepath)


def anndata_pbmc3k_processed() -> ad.AnnData:  # pragma: no cover
    """Modified from scanpy.pbmc3k_processed()."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/scrnaseq_scanpy_pbmc3k_processed.h5ad"
    )
    pbmc3k = ad.read_h5ad(filepath)
    pbmc3k.obs.rename(columns={"louvain": "cell_type"}, inplace=True)
    return pbmc3k


def anndata_human_immune_cells(
    populate_registries: bool = False,
) -> ad.AnnData:  # pragma: no cover
    """Cross-tissue immune cell analysis reveals tissue-specific features in humans.

    From: https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3  # noqa
    Dataset: Global

    To reproduce the subsample::

        adata = sc.read('Global.h5ad')
        adata.obs = adata.obs[['donor_id', 'tissue', 'cell_type', 'assay', 'tissue_ontology_term_id', 'cell_type_ontology_term_id', 'assay_ontology_term_id']].copy()
        sc.pp.subsample(adata, fraction=0.005)
        del adata.uns["development_stage_ontology_term_id_colors"]
        del adata.uns["sex_ontology_term_id_colors"]
        adata.write('human_immune.h5ad')
    """
    filepath, _ = urlretrieve("https://lamindb-test.s3.amazonaws.com/human_immune.h5ad")
    adata = ad.read_h5ad(filepath)
    adata.var.drop(columns=["gene_symbols", "feature_name"], inplace=True)
    adata.uns.pop("cell_type_ontology_term_id_colors")
    adata.uns.pop("title")
    adata.uns.pop("schema_version")
    adata.obs.columns = adata.obs.columns.str.replace("donor_id", "donor")
    columns = [col for col in adata.obs.columns if "ontology_term" not in col]
    adata.obs = adata.obs[columns]
    if populate_registries:
        import lnschema_bionty as lb

        import lamindb as ln

        verbosity = ln.settings.verbosity
        ln.settings.verbosity = "error"
        auto_save_parents = lb.settings.auto_save_parents
        lb.settings.auto_save_parents = False
        ln.save(
            lb.Gene.from_values(
                adata.var.index, field="ensembl_gene_id", organism="human"
            )
        )
        ln.save(lb.CellType.from_values(adata.obs.cell_type, field="name"))
        ln.save(lb.ExperimentalFactor.from_values(adata.obs.assay, field="name"))
        ln.save(lb.Tissue.from_values(adata.obs.tissue, field="name"))
        ln.Feature(name="cell_type", type="category", registries=[lb.CellType]).save()
        ln.Feature(
            name="assay", type="category", registries=[lb.ExperimentalFactor]
        ).save()
        ln.Feature(name="tissue", type="category", registries=[lb.Tissue]).save()
        ln.Feature(name="organism", type="category", registries=[lb.Organism]).save()
        ln.Feature(name="donor", type="category", registries=[ln.ULabel]).save()
        lb.ExperimentalFactor.from_bionty(ontology_id="EFO:0008913").save()
        ln.save([ln.ULabel(name=name) for name in adata.obs.donor.unique()])
        ln.settings.verbosity = verbosity
        lb.settings.auto_save_parents = auto_save_parents
    return adata


def anndata_with_obs() -> ad.AnnData:
    """Create a mini anndata with cell_type, disease and tissue."""
    import anndata as ad
    import bionty as bt

    celltypes = ["T cell", "hematopoietic stem cell", "hepatocyte", "my new cell type"]
    celltype_ids = ["CL:0000084", "CL:0000037", "CL:0000182", ""]
    diseases = [
        "chronic kidney disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
        "Alzheimer disease",
    ]
    tissues = ["kidney", "liver", "heart", "brain"]
    df = pd.DataFrame()
    df["cell_type"] = celltypes * 10
    df["cell_type_id"] = celltype_ids * 10
    df["tissue"] = tissues * 10
    df["disease"] = diseases * 10
    df.index = "obs" + df.index.astype(str)

    adata = ad.AnnData(X=np.zeros(shape=(40, 100), dtype=np.float32), obs=df)
    adata.var.index = bt.Gene().df().head(100)["ensembl_gene_id"].values

    return adata


def anndata_suo22_Visium10X():  # pragma: no cover
    """AnnData from Suo22 generated by 10x Visium."""
    import anndata as ad

    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/suo22_Visium10X_data_LI_subset.h5ad",
        "Visium10X_data_LI_subset.h5ad",
    )
    Path("suo22/").mkdir(exist_ok=True)
    filepath = Path(filepath).rename("suo22/Visium10X_data_LI_subset.h5ad")
    return ad.read_h5ad(filepath)


def mudata_papalexi21_subset():  # pragma: no cover
    """A subsetted mudata from papalexi21.

    To reproduce the subsetting:
    >>> !wget https://figshare.com/ndownloader/files/36509460
    >>> import mudata as md
    >>> import scanpy as sc
    >>> mdata = md.read_h5mu("36509460")
    >>> mdata = sc.pp.subsample(mdata, n_obs=200, copy=True)[0]
    >>> mdata[:, -300:].copy().write("papalexi21_subset_200x300_lamindb_demo_2023-07-25.h5mu")  # noqa
    """
    import mudata as md

    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/papalexi21_subset_200x300_lamindb_demo_2023-07-25.h5mu",  # noqa
        "papalexi21_subset.h5mu",
    )

    return md.read_h5mu(filepath)


def df_iris() -> pd.DataFrame:
    """The iris dataset as in sklearn.

    Original code::

        sklearn.datasets.load_iris(as_frame=True).frame
    """
    filepath, _ = urlretrieve("https://lamindb-test.s3.amazonaws.com/iris.parquet")
    return pd.read_parquet(filepath)


def df_iris_in_meter() -> pd.DataFrame:
    """The iris dataset with lengths in meter."""
    df = df_iris()
    # rename columns
    df.rename(
        columns={
            "sepal length (cm)": "sepal_length",
            "sepal width (cm)": "sepal_width",
            "petal length (cm)": "petal_length",
            "petal width (cm)": "petal_width",
        },
        inplace=True,
    )
    df[["sepal_length", "sepal_width", "petal_length", "petal_width"]] /= 100
    df["iris_organism_name"] = df["target"].map(
        {0: "setosa", 1: "versicolor", 2: "virginica"}
    )
    del df["target"]
    return df


def df_iris_in_meter_study1() -> pd.DataFrame:
    """The iris dataset with lengths in meter."""
    df_iris = df_iris_in_meter()
    return df_iris.iloc[: len(df_iris) // 2]


def df_iris_in_meter_study2() -> pd.DataFrame:
    """The iris dataset with lengths in meter."""
    df_iris = df_iris_in_meter()
    return df_iris.iloc[len(df_iris) // 2 :]


def dir_scrnaseq_cellranger(
    sample_name: str, basedir: Union[str, Path] = "./", output_only: bool = True
):  # pragma: no cover
    """Generate mock cell ranger outputs.

    Args:
        sample_name: name of the sample
        basedir: run directory
        output_only: only generate output files
    """
    basedir = Path(basedir)

    if not output_only:  # pragma: no cover
        fastqdir = basedir / "fastq"
        fastqdir.mkdir(parents=True, exist_ok=True)
        fastqfile1 = fastqdir / f"{sample_name}_R1_001.fastq.gz"
        with open(fastqfile1, "w") as f:
            f.write(f"{ids.base62(n_char=6)}")
        fastqfile2 = fastqdir / f"{sample_name}_R2_001.fastq.gz"
        fastqfile2.touch(exist_ok=True)
        with open(fastqfile2, "w") as f:
            f.write(f"{ids.base62(n_char=6)}")

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
            f.write(f"{ids.base62(n_char=6)}")

    return sampledir


def schmidt22_crispra_gws_IFNG(basedir=".") -> Path:  # pragma: no cover
    """CRISPRi screen dataset of Schmidt22.

    Originally from: https://zenodo.org/record/5784651
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/schmidt22-crispra-gws-IFNG.csv",
        "schmidt22-crispra-gws-IFNG.csv",
    )
    return Path(filepath).rename(Path(basedir) / filepath)


def schmidt22_perturbseq(basedir=".") -> Path:  # pragma: no cover
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
    return Path(filepath).rename(Path(basedir) / filepath)
