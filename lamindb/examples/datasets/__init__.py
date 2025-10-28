"""Example datasets.

The mini immuno dataset
-----------------------

.. autofunction:: mini_immuno

Small in-memory datasets
------------------------

.. autofunction:: anndata_with_obs

Files
-----

.. autofunction:: file_fcs
.. autofunction:: file_fcs_alpert19
.. autofunction:: file_tsv_rnaseq_nfcore_salmon_merged_gene_counts
.. autofunction:: file_jpg_paradisi05
.. autofunction:: file_tiff_suo22
.. autofunction:: file_fastq
.. autofunction:: file_bam
.. autofunction:: file_mini_csv

Directories
-----------

.. autofunction:: dir_scrnaseq_cellranger
.. autofunction:: dir_iris_images

Dictionary, Dataframe, AnnData, MuData, SpatialData
----------------------------------------------------

.. autofunction:: dict_cellxgene_uns
.. autofunction:: df_iris
.. autofunction:: df_iris_in_meter
.. autofunction:: df_iris_in_meter_study1
.. autofunction:: df_iris_in_meter_study2
.. autofunction:: anndata_mouse_sc_lymph_node
.. autofunction:: anndata_human_immune_cells
.. autofunction:: anndata_pbmc68k_reduced
.. autofunction:: anndata_file_pbmc68k_test
.. autofunction:: anndata_pbmc3k_processed
.. autofunction:: anndata_with_obs
.. autofunction:: anndata_suo22_Visium10X
.. autofunction:: mudata_papalexi21_subset
.. autofunction:: schmidt22_crispra_gws_IFNG
.. autofunction:: schmidt22_perturbseq
.. autofunction:: spatialdata_blobs

Other
-----

.. autofunction:: fake_bio_notebook_titles
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
    dict_cellxgene_uns,
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
