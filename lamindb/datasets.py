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
