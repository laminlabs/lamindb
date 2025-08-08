"""Example Croissant files.

Examples for MLCommons Croissant files, which are used to store metadata about datasets.
"""

import json
from pathlib import Path


def mini_immuno(n_files: int = 1) -> list[Path]:
    """Return paths to the mini immuno dataset and its metadata as a Croissant file.

    Args:
        n_files: Number of files inside the croissant file. Default is 1.
    """
    from ..datasets import file_mini_csv
    from ..datasets.mini_immuno import get_dataset1

    adata = get_dataset1(otype="AnnData")
    dataset1_path = Path("mini_immuno.anndata.zarr")
    adata.write_zarr(dataset1_path)
    orig_croissant_path = (
        Path(__file__).parent / "mini_immuno.anndata.zarr_metadata.json"
    )
    with open(orig_croissant_path, encoding="utf-8") as f:
        data = json.load(f)
    if n_files == 2:
        dataset2_path = file_mini_csv()
        data["distribution"].append(
            {
                "@type": "sc:FileObject",
                "@id": "mini.csv",
                "name": "mini.csv",
                "encodingFormat": "text/csv",
            }
        )
    croissant_path = Path("mini_immuno.anndata.zarr_metadata.json")
    with open(croissant_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    result: list[Path] = [croissant_path, dataset1_path]
    if n_files == 1:
        return result
    result.append(dataset2_path)
    return result
