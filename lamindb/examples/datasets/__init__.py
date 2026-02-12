"""Example datasets.

The mini immuno dataset
-----------------------

.. autosummary::
   :toctree: .

   mini_immuno

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
.. autofunction:: anndata_suo22_Visium10X
.. autofunction:: anndata_visium_mouse_cellxgene
.. autofunction:: mudata_papalexi21_subset
.. autofunction:: schmidt22_crispra_gws_IFNG
.. autofunction:: schmidt22_perturbseq
.. autofunction:: spatialdata_blobs


Other
-----

.. autofunction:: fake_bio_notebook_titles
"""

import importlib.util
import sys


def __getattr__(name: str):
    """Lazy-import datasets to avoid loading pandas/anndata at package import."""
    if name == "mini_immuno":
        # Use importlib to avoid __getattr__ recursion when importing submodule
        spec = importlib.util.find_spec(
            "lamindb.examples.datasets.mini_immuno",
            package="lamindb.examples.datasets",
        )
        if spec is None or spec.loader is None:
            raise ImportError("Could not find module mini_immuno")
        module = importlib.util.module_from_spec(spec)
        sys.modules["lamindb.examples.datasets.mini_immuno"] = module
        spec.loader.exec_module(module)
        return module
    if name in ("small_dataset1", "small_dataset2"):
        mini_immuno = importlib.import_module(
            ".mini_immuno", package="lamindb.examples.datasets"
        )
        return (
            mini_immuno.get_dataset1
            if name == "small_dataset1"
            else mini_immuno.get_dataset2
        )
    _core_names = (
        "anndata_file_pbmc68k_test",
        "anndata_human_immune_cells",
        "anndata_mouse_sc_lymph_node",
        "anndata_pbmc3k_processed",
        "anndata_pbmc68k_reduced",
        "anndata_suo22_Visium10X",
        "df_iris",
        "df_iris_in_meter",
        "df_iris_in_meter_study1",
        "df_iris_in_meter_study2",
        "dict_cellxgene_uns",
        "dir_iris_images",
        "dir_scrnaseq_cellranger",
        "file_bam",
        "file_fastq",
        "file_fcs",
        "file_fcs_alpert19",
        "file_jpg_paradisi05",
        "file_mini_csv",
        "file_tiff_suo22",
        "file_tsv_rnaseq_nfcore_salmon_merged_gene_counts",
        "mudata_papalexi21_subset",
        "schmidt22_crispra_gws_IFNG",
        "schmidt22_perturbseq",
        "spatialdata_blobs",
        "anndata_visium_mouse_cellxgene",
    )
    if name in _core_names:
        _core = importlib.import_module("._core", package="lamindb.examples.datasets")
        return getattr(_core, name)
    if name in ("anndata_with_obs", "small_dataset3_cellxgene"):
        _small = importlib.import_module("._small", package="lamindb.examples.datasets")
        return getattr(_small, name)
    if name == "fake_bio_notebook_titles":
        _fake = importlib.import_module("._fake", package="lamindb.examples.datasets")
        return _fake.fake_bio_notebook_titles
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
