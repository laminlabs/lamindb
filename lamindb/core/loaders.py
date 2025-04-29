"""Loaders in :class:`lamindb.Artifact.load`.

.. autosummary::
   :toctree: .

   SUPPORTED_SUFFIXES
   load_fcs
   load_tsv
   load_h5ad
   load_h5mu
   load_html
   load_json
   load_image
   load_svg

"""

from __future__ import annotations

import builtins
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from anndata import read_h5ad
from lamin_utils import logger
from lamindb_setup.core.upath import (
    create_path,
    infer_filesystem,
)

from ..core._settings import settings

if TYPE_CHECKING:
    from anndata import AnnData
    from lamindb_setup.core.types import UPathStr
    from mudata import MuData

    from lamindb.core.types import ScverseDataStructures

try:
    from ..core.storage._zarr import load_zarr
except ImportError:

    def load_zarr(storepath):  # type: ignore
        raise ImportError("Please install zarr: pip install 'lamindb[zarr]'")


is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


# tested in lamin-usecases
def load_fcs(*args, **kwargs) -> AnnData:
    """Load an `.fcs` file to `AnnData`."""
    try:
        import readfcs
    except ImportError:  # pragma: no cover
        raise ImportError("Please install readfcs: pip install readfcs") from None
    return readfcs.read(*args, **kwargs)


def load_tsv(path: UPathStr, **kwargs) -> pd.DataFrame:
    """Load `.tsv` file to `DataFrame`."""
    path_sanitized = Path(path)
    return pd.read_csv(path_sanitized, sep="\t", **kwargs)


def load_h5ad(filepath, **kwargs) -> AnnData:
    """Load an `.h5ad` file to `AnnData`."""
    fs, filepath = infer_filesystem(filepath)
    compression = kwargs.pop("compression", "infer")
    with fs.open(filepath, mode="rb", compression=compression) as file:
        adata = read_h5ad(file, backed=False, **kwargs)
        return adata


def load_h5mu(filepath: UPathStr, **kwargs) -> MuData:
    """Load an `.h5mu` file to `MuData`."""
    import mudata as md

    path_sanitized = Path(filepath)
    return md.read_h5mu(path_sanitized, **kwargs)


def load_html(path: UPathStr) -> None | UPathStr:
    """Display `.html` in ipython, otherwise return path."""
    if is_run_from_ipython:
        with open(path, encoding="utf-8") as f:
            html_content = f.read()
        # Extract the body content using regular expressions
        body_content = re.findall(
            r"<body(?:.*?)>(?:.*?)</body>", html_content, re.DOTALL
        )
        # Remove any empty body tags
        if body_content:
            body_content = body_content[0]
            body_content = body_content.strip()  # type: ignore
        from IPython.display import HTML, display

        display(HTML(data=body_content))
        return None
    else:
        return path


def load_json(path: UPathStr) -> dict[str, Any] | list[Any]:
    """Load `.json` to `dict`."""
    import json

    with open(path) as f:
        data = json.load(f)
    return data


def load_yaml(path: UPathStr) -> dict[str, Any] | list[Any]:
    """Load `.yaml` to `dict`."""
    import yaml  # type: ignore

    with open(path) as f:
        data = yaml.safe_load(f)
    return data


def load_image(path: UPathStr) -> None | UPathStr:
    """Display `.jpg`, `.gif` or `.png` in ipython, otherwise return path."""
    if is_run_from_ipython:
        from IPython.display import Image, display

        display(Image(filename=path))
        return None
    else:
        return path


def load_svg(path: UPathStr) -> None | UPathStr:
    """Display `.svg` in ipython, otherwise return path."""
    if is_run_from_ipython:
        from IPython.display import SVG, display

        display(SVG(filename=path))
        return None
    else:
        return path


def load_rds(path: UPathStr) -> UPathStr:
    """Just warn when trying to load `.rds`."""
    logger.warning("Please use `laminr` to load `.rds` files")
    return path


FILE_LOADERS = {
    ".csv": pd.read_csv,
    ".csv.gz": pd.read_csv,
    ".tsv": load_tsv,
    ".tsv.gz": load_tsv,
    ".h5ad": load_h5ad,
    ".h5ad.gz": load_h5ad,
    ".parquet": pd.read_parquet,
    ".parquet.gz": pd.read_parquet,  # this doesn't work for externally gzipped files, REMOVE LATER
    ".fcs": load_fcs,
    ".zarr": load_zarr,
    ".html": load_html,
    ".json": load_json,
    ".yaml": load_yaml,
    ".h5mu": load_h5mu,
    ".gif": load_image,
    ".jpg": load_image,
    ".png": load_image,
    ".svg": load_svg,
    ".rds": load_rds,
}

SUPPORTED_SUFFIXES = [sfx for sfx in FILE_LOADERS.keys() if sfx != ".rds"]
"""Suffixes with defined artifact loaders."""


def load_to_memory(
    filepath: UPathStr, **kwargs
) -> (
    pd.DataFrame | ScverseDataStructures | dict[str, Any] | list[Any] | UPathStr | None
):
    """Load a file into memory.

    Returns the filepath if no in-memory form is found.
    May return None in interactive sessions for images.
    """
    filepath = create_path(filepath)

    # infer the correct suffix when .gz is present
    suffixes = filepath.suffixes
    suffix = (
        "".join(suffixes[-2:])
        if len(suffixes) > 1 and ".gz" in suffixes
        else filepath.suffix
    )

    loader = FILE_LOADERS.get(suffix, None)
    if loader is None:
        raise NotImplementedError(
            f"There is no loader for {suffix} files. Use .cache() to get the path."
        )

    filepath = settings._storage_settings.cloud_to_local(filepath, print_progress=True)

    return loader(filepath, **kwargs)
