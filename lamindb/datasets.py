"""Small example datasets."""

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


def bfx_output() -> Path:
    """Directory with exemplary BFX output."""
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/bfx-output.zip",
    )
    from zipfile import ZipFile

    with ZipFile(filepath, "r") as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall(path=".")

    return Path("bfx-output")


def file_mouse_sc_lymph_node() -> Path:
    """Mouse lymph node scRNA-seq dataset from EBI.

    Subsampled to 10k genes.

    From: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8414/
    """
    filepath, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/E-MTAB-8414.h5ad",  # noqa
        "mouse_sc_lymph_node.h5ad",
    )
    return Path(filepath)
