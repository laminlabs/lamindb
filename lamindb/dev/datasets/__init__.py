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
   anndata_mouse_sc_lymph_node
   anndata_human_immune_cells
   anndata_pbmc68k_reduced
   anndata_pbmc3k_processed
"""

from ._core import (
    anndata_human_immune_cells,
    anndata_mouse_sc_lymph_node,
    anndata_pbmc3k_processed,
    anndata_pbmc68k_reduced,
    dir_scrnaseq_cellranger,
    file_bam,
    file_fastq,
    file_fcs,
    file_jpg_paradisi05,
    file_mini_csv,
    generate_cell_ranger_files,
)
