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

import builtins
import re
from pathlib import Path
from typing import TYPE_CHECKING

import anndata as ad
import pandas as pd
from lamindb_setup.core.types import UPathStr
from lamindb_setup.core.upath import (
    create_path,
    infer_filesystem,
)

from lamindb.core._settings import settings

if TYPE_CHECKING:
    import mudata as md

try:
    from .storage._zarr import load_anndata_zarr
except ImportError:

    def load_anndata_zarr(storepath):  # type: ignore
        raise ImportError("Please install zarr: pip install zarr")


is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


# tested in lamin-usecases
def load_fcs(*args, **kwargs) -> ad.AnnData:
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


def load_h5ad(filepath, **kwargs) -> ad.AnnData:
    """Load an `.h5ad` file to `AnnData`."""
    fs, filepath = infer_filesystem(filepath)

    with fs.open(filepath, mode="rb") as file:
        adata = ad.read_h5ad(file, backed=False, **kwargs)
        return adata


def load_h5mu(filepath: UPathStr, **kwargs):
    """Load an `.h5mu` file to `MuData`."""
    import mudata as md

    path_sanitized = Path(filepath)
    return md.read_h5mu(path_sanitized, **kwargs)


def load_html(path: UPathStr):
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
    else:
        return path


def load_json(path: UPathStr) -> dict:
    """Load `.json` to `dict`."""
    import json

    with open(path) as f:
        data = json.load(f)
    return data


def load_image(path: UPathStr):
    """Display `.svg` in ipython, otherwise return path."""
    if is_run_from_ipython:
        from IPython.display import Image, display

        display(Image(filename=path))
    else:
        return path


def load_svg(path: UPathStr) -> None | Path:
    """Display `.svg` in ipython, otherwise return path."""
    if is_run_from_ipython:
        from IPython.display import SVG, display

        display(SVG(filename=path))
        return None
    else:
        return path


FILE_LOADERS = {
    ".csv": pd.read_csv,
    ".tsv": load_tsv,
    ".h5ad": load_h5ad,
    ".parquet": pd.read_parquet,
    ".fcs": load_fcs,
    ".zarr": load_anndata_zarr,
    ".html": load_html,
    ".json": load_json,
    ".h5mu": load_h5mu,
    ".jpg": load_image,
    ".png": load_image,
    ".svg": load_svg,
}

SUPPORTED_SUFFIXES = list(FILE_LOADERS.keys())
"""Suffixes with defined artifact loaders."""


def load_to_memory(filepath: UPathStr, stream: bool = False, **kwargs):
    """Load a file into memory.

    Returns the filepath if no in-memory form is found.
    """
    filepath = create_path(filepath)

    if filepath.suffix not in {".h5ad", ".zarr"}:
        stream = False

    if not stream:
        # caching happens here if filename is a UPath
        # todo: make it safe when filepath is just Path
        filepath = settings._storage_settings.cloud_to_local(
            filepath, print_progress=True
        )

    loader = FILE_LOADERS.get(filepath.suffix)
    if loader is None:
        return filepath
    else:
        return loader(filepath, **kwargs)
