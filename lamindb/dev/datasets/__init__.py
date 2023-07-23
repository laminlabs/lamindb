"""Small example datasets.

.. autosummary::
   :toctree: .

   file_fcs
   file_jpg_paradisi05
   file_fastq
   file_bam
   file_mini_csv
   dir_scrnaseq_cellranger
   generate_cell_ranger_files
   df_iris
   df_iris_in_meter
   df_iris_in_meter_batch1
   df_iris_in_meter_batch2
   anndata_mouse_sc_lymph_node
   anndata_human_immune_cells
   anndata_pbmc68k_reduced
   anndata_pbmc3k_processed
   anndata_with_obs
   schmidt22_crispra_gws_IFNG
   schmidt22_perturbseq
   fake_bio_notebook_titles
"""

from ._core import (
    anndata_human_immune_cells,
    anndata_mouse_sc_lymph_node,
    anndata_pbmc3k_processed,
    anndata_pbmc68k_reduced,
    anndata_with_obs,
    df_iris,
    df_iris_in_meter,
    df_iris_in_meter_batch1,
    df_iris_in_meter_batch2,
    dir_scrnaseq_cellranger,
    file_bam,
    file_fastq,
    file_fcs,
    file_jpg_paradisi05,
    file_mini_csv,
    generate_cell_ranger_files,
    schmidt22_crispra_gws_IFNG,
    schmidt22_perturbseq,
)
from ._fake import fake_bio_notebook_titles
