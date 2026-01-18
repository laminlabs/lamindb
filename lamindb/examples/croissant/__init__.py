"""Examples for MLCommons Croissant files, which are used to store metadata about datasets.

.. autofunction:: mini_immuno

"""

import json
from pathlib import Path


def mini_immuno(
    n_files: int = 1, filepath_prefix: str = "", strip_version: bool = False
) -> list[Path]:
    """Return paths to the mini immuno dataset and its metadata as a Croissant file.

    Args:
        n_files: Number of files inside the croissant file.
        filepath_prefix: Move the dataset and references to it in a specific directory.

    Example

        ::

            croissant_path, dataset1_path = ln.examples.croissant.mini_immuno()
            croissant_path, dataset1_path, dataset2_path = ln.examples.croissant.mini_immuno(n_files=2)
    """
    from ..datasets import file_mini_csv
    from ..datasets.mini_immuno import get_dataset1

    adata = get_dataset1(otype="AnnData")
    if filepath_prefix:
        dataset1_path = Path(filepath_prefix) / "mini_immuno.anndata.zarr"
    else:
        dataset1_path = Path("mini_immuno.anndata.zarr")
    adata.write_zarr(dataset1_path)
    orig_croissant_path = (
        Path(__file__).parent / "mini_immuno.anndata.zarr_metadata.json"
    )
    with open(orig_croissant_path, encoding="utf-8") as f:
        data = json.load(f)
    if filepath_prefix:
        assert data["distribution"][0]["@id"] == "mini_immuno.anndata.zarr"  # noqa: S101
        data["distribution"][0]["@id"] = str(Path(filepath_prefix) / dataset1_path.name)
    if strip_version:
        data.pop("version", None)
    if n_files == 2:
        file_mini_csv()
        if filepath_prefix:
            dataset2_path = Path(filepath_prefix) / "mini.csv"
        else:
            dataset2_path = Path("mini.csv")
        data["distribution"].append(
            {
                "@type": "sc:FileObject",
                "@id": dataset2_path.as_posix(),
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
