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


def folder_bfx_output() -> Path:
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/bfx-output.zip",
    )
    from zipfile import ZipFile

    with ZipFile(filepath, "r") as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall(path=".")

    return Path("bfx-output")


def file_10x_mouse_sc_E18_brain() -> Path:
    """10x dataset - 1k Brain Cells from an E18 Mouse (v3 chemistry).

    Subsampled to 5k genes.

    From: https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0  # noqa
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/neuron_1k_v3_filtered_feature_bc_matrix.h5ad",  # noqa
        "mouse_sc_E18_brain.h5ad",
    )
    return Path(filepath)
