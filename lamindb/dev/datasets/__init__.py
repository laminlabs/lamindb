"""Test datasets.

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
   dir_scrnaseq_cellranger
   dir_iris_images
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
   fake_bio_notebook_titles
"""

from ._core import (
    anndata_file_pbmc68k_test,
    anndata_human_immune_cells,
    anndata_mouse_sc_lymph_node,
    anndata_pbmc3k_processed,
    anndata_pbmc68k_reduced,
    anndata_suo22_Visium10X,
    anndata_with_obs,
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
)
from ._fake import fake_bio_notebook_titles
