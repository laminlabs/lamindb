"""Test datasets.

The mini immuno dataset.

.. autosummary::
   :toctree: .

   mini_immuno

Small in-memory datasets.

.. autosummary::
   :toctree: .

   anndata_with_obs

Files.

.. autosummary::
   :toctree: .

   file_fcs
   file_fcs_alpert19
   file_tsv_rnaseq_nfcore_salmon_merged_gene_counts
   file_jpg_paradisi05
   file_tiff_suo22
   file_fastq
   file_bam
   file_mini_csv

Directories.

.. autosummary::
   :toctree: .

   dir_scrnaseq_cellranger
   dir_iris_images

Dataframe, AnnData, MuData.

.. autosummary::
   :toctree: .

   df_iris
   df_iris_in_meter
   df_iris_in_meter_study1
   df_iris_in_meter_study2
   anndata_mouse_sc_lymph_node
   anndata_human_immune_cells
   anndata_pbmc68k_reduced
   anndata_file_pbmc68k_test
   anndata_pbmc3k_processed
   anndata_with_obs
   anndata_suo22_Visium10X
   mudata_papalexi21_subset
   schmidt22_crispra_gws_IFNG
   schmidt22_perturbseq

Other.

.. autosummary::
   :toctree: .

   fake_bio_notebook_titles
"""

from . import mini_immuno
from ._core import (
    anndata_file_pbmc68k_test,
    anndata_human_immune_cells,
    anndata_mouse_sc_lymph_node,
    anndata_pbmc3k_processed,
    anndata_pbmc68k_reduced,
    anndata_suo22_Visium10X,
    df_iris,
    df_iris_in_meter,
    df_iris_in_meter_study1,
    df_iris_in_meter_study2,
    dir_iris_images,
    dir_scrnaseq_cellranger,
    file_bam,
    file_fastq,
    file_fcs,
    file_fcs_alpert19,
    file_jpg_paradisi05,
    file_mini_csv,
    file_tiff_suo22,
    file_tsv_rnaseq_nfcore_salmon_merged_gene_counts,
    mudata_papalexi21_subset,
    schmidt22_crispra_gws_IFNG,
    schmidt22_perturbseq,
    spatialdata_blobs,
)
from ._fake import fake_bio_notebook_titles
from ._small import (
    anndata_with_obs,
    small_dataset3_cellxgene,
)

small_dataset1 = mini_immuno.get_dataset1  # backward compat
small_dataset2 = mini_immuno.get_dataset2  # backward compat
