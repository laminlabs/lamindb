from pathlib import Path
from urllib.request import urlretrieve


def fcs_file() -> Path:
    """Return fcs file example."""
    fcs_filename, _ = urlretrieve(
        "https://lamindb-test.s3.amazonaws.com/example.fcs", "example.fcs"
    )
    return Path(fcs_filename)
