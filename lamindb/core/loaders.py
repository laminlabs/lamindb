"""Loaders in :class:`lamindb.Artifact.load`.

.. autodata:: SUPPORTED_SUFFIXES
.. autofunction:: load_fcs
.. autofunction:: load_tsv
.. autofunction:: load_h5ad
.. autofunction:: load_h5mu
.. autofunction:: load_html
.. autofunction:: load_json
.. autofunction:: load_image
.. autofunction:: load_svg

"""

from __future__ import annotations

import builtins
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core.upath import (
    create_path,
    extract_suffix_from_path,
    infer_filesystem,
)

if TYPE_CHECKING:
    import pandas as pd
    from anndata import AnnData
    from lamindb_setup.types import UPathStr
    from mudata import MuData

    from lamindb.core.types import ScverseDataStructures


def load_zarr(storepath, **kwargs):
    """Lazy-import to avoid loading storage at package import."""
    try:
        from ..core.storage._zarr import load_zarr as _load_zarr

        return _load_zarr(storepath, **kwargs)
    except ImportError:
        raise ImportError("Please install zarr: pip install 'lamindb[zarr]'") from None


is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


# tested in lamin-usecases
def load_fcs(*args, **kwargs) -> AnnData:
    """Load an `.fcs` file to `AnnData`."""
    try:
        import readfcs
    except ImportError:  # pragma: no cover
        raise ImportError("Please install readfcs: pip install readfcs") from None
    return readfcs.read(*args, **kwargs)


def load_tsv(path: UPathStr, **kwargs):
    """Load `.tsv` file to `DataFrame`."""
    import pandas as pd

    path_sanitized = Path(path)
    return pd.read_csv(path_sanitized, sep="\t", **kwargs)


def load_h5ad(filepath, **kwargs) -> AnnData:
    """Load an `.h5ad` file to `AnnData`."""
    from anndata import read_h5ad

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


def load_txt(path: Path) -> str:
    """Load `.txt` file to `str`."""
    return path.read_text(encoding="utf-8")


def load_rds(path: UPathStr) -> UPathStr:
    """Just warn when trying to load `.rds`."""
    logger.warning("Please use `laminr` to load `.rds` files")
    return path


_file_loaders_cache: dict[str, Any] | None = None


def _get_file_loaders() -> dict[str, Any]:
    """Lazy-build FILE_LOADERS to avoid importing pandas/anndata at module load."""
    global _file_loaders_cache
    if _file_loaders_cache is None:
        import pandas as pd

        _file_loaders_cache = {
            ".csv": pd.read_csv,
            ".csv.gz": pd.read_csv,
            ".csv.tar.gz": pd.read_csv,
            ".tsv": load_tsv,
            ".tsv.gz": load_tsv,
            ".tsv.tar.gz": load_tsv,
            ".h5ad": load_h5ad,
            ".h5ad.gz": load_h5ad,
            ".h5ad.tar.gz": load_h5ad,
            ".parquet": pd.read_parquet,
            ".fcs": load_fcs,
            ".zarr": load_zarr,
            ".anndata.zarr": load_zarr,
            ".html": load_html,
            ".json": load_json,
            ".vitessce.json": load_json,
            ".yaml": load_yaml,
            ".h5mu": load_h5mu,
            ".gif": load_image,
            ".jpg": load_image,
            ".png": load_image,
            ".svg": load_svg,
            ".rds": load_rds,
            ".txt": load_txt,
            ".fasta": load_txt,
        }
    return _file_loaders_cache


SUPPORTED_SUFFIXES = [
    ".csv",
    ".csv.gz",
    ".csv.tar.gz",
    ".tsv",
    ".tsv.gz",
    ".tsv.tar.gz",
    ".h5ad",
    ".h5ad.gz",
    ".h5ad.tar.gz",
    ".parquet",
    ".fcs",
    ".zarr",
    ".anndata.zarr",
    ".html",
    ".json",
    ".vitessce.json",
    ".yaml",
    ".h5mu",
    ".gif",
    ".jpg",
    ".png",
    ".svg",
    ".txt",
    ".fasta",
]
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
    suffix = extract_suffix_from_path(filepath)
    loader = _get_file_loaders().get(suffix, None)
    if loader is None:
        raise NotImplementedError(
            f"There is no loader for {suffix} files. Use .cache() to get the path."
        )

    filepath = setup_settings.paths.cloud_to_local(filepath, print_progress=True)

    return loader(filepath, **kwargs)
