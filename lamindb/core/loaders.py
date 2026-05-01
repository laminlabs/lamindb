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
from typing import TYPE_CHECKING, Any, Callable, cast

from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core.upath import (
    create_path,
    extract_suffix_from_path,
    infer_filesystem,
)

if TYPE_CHECKING:
    from anndata import AnnData
    from lamindb_setup.types import AnyPathStr
    from mudata import MuData
    from pandas import DataFrame

    from lamindb.core.storage.types import ScverseDataStructures


is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


# tested in lamin-usecases
def load_fcs(*args, **kwargs) -> AnnData:
    """Load an `.fcs` file to `AnnData`."""
    try:
        import readfcs
    except ImportError:  # pragma: no cover
        raise ImportError("Please install readfcs: pip install readfcs") from None
    return readfcs.read(*args, **kwargs)


# for types below note that local UPaths are subclasses of Path
# so Path(UPath(...)) properly coerces local UPaths and throws an error for cloud UPaths


def load_csv(path: AnyPathStr, **kwargs) -> DataFrame:
    """Load `.csv` file to `DataFrame`."""
    import pandas as pd

    path_sanitized = Path(path)
    return pd.read_csv(path_sanitized, **kwargs)


def load_parquet(path: AnyPathStr, **kwargs) -> DataFrame:
    """Load `.parquet` file to `DataFrame`."""
    import pandas as pd

    path_sanitized = Path(path)
    return pd.read_parquet(path_sanitized, **kwargs)


def load_tsv(path: AnyPathStr, **kwargs) -> DataFrame:
    """Load `.tsv` file to `DataFrame`."""
    import pandas as pd

    path_sanitized = Path(path)
    return pd.read_csv(path_sanitized, sep="\t", **kwargs)


def load_h5ad(filepath: AnyPathStr, **kwargs) -> AnnData:
    """Load an `.h5ad` file to `AnnData`."""
    from anndata import read_h5ad

    fs, filepath_str = infer_filesystem(filepath)
    compression = kwargs.pop("compression", "infer")
    with fs.open(filepath_str, mode="rb", compression=compression) as file:
        adata = read_h5ad(file, backed=False, **kwargs)
        return adata


def load_h5mu(filepath: AnyPathStr, **kwargs) -> MuData:
    """Load an `.h5mu` file to `MuData`."""
    import mudata as md

    path_sanitized = Path(filepath)
    return md.read_h5mu(path_sanitized, **kwargs)


def load_zarr(storepath, **kwargs):  # type: ignore
    try:
        from ..core.storage._zarr import load_zarr as _load_zarr
    except ImportError:
        raise ImportError("Please install zarr: pip install 'lamindb[zarr]'") from None
    return _load_zarr(storepath, **kwargs)


def load_html(path: AnyPathStr) -> None | AnyPathStr:
    """Display `.html` in ipython, otherwise return path."""
    if is_run_from_ipython:
        path_sanitized = Path(path)
        with path_sanitized.open(encoding="utf-8") as f:
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


def load_json(path: AnyPathStr) -> dict[str, Any] | list[Any]:
    """Load `.json` to `dict`."""
    import json

    path_sanitized = Path(path)
    with path_sanitized.open(encoding="utf-8") as f:
        data = json.load(f)
    return data


def load_yaml(path: AnyPathStr) -> dict[str, Any] | list[Any]:
    """Load `.yaml` to `dict`."""
    import yaml  # type: ignore

    path_sanitized = Path(path)
    with path_sanitized.open(encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data


def load_image(path: AnyPathStr) -> None | AnyPathStr:
    """Display `.jpg`, `.gif` or `.png` in ipython, otherwise return path."""
    if is_run_from_ipython:
        from IPython.display import Image, display

        path_sanitized = Path(path)
        display(Image(filename=path_sanitized.as_posix()))
        return None
    else:
        return path


def load_svg(path: AnyPathStr) -> None | AnyPathStr:
    """Display `.svg` in ipython, otherwise return path."""
    if is_run_from_ipython:
        from IPython.display import SVG, display

        path_sanitized = Path(path)
        display(SVG(filename=path_sanitized.as_posix()))
        return None
    else:
        return path


def load_txt(path: AnyPathStr) -> str:
    """Load `.txt` file to `str`."""
    path_sanitized = Path(path)
    return path_sanitized.read_text(encoding="utf-8")


def load_rds(path: AnyPathStr) -> AnyPathStr:
    """Just warn when trying to load `.rds`."""
    logger.warning("Please use `laminr` to load `.rds` files")
    return path


FILE_LOADERS = {
    ".csv": load_csv,
    ".csv.gz": load_csv,
    ".csv.tar.gz": load_csv,
    ".tsv": load_tsv,
    ".tsv.gz": load_tsv,
    ".tsv.tar.gz": load_tsv,
    ".h5ad": load_h5ad,
    ".h5ad.gz": load_h5ad,
    ".h5ad.tar.gz": load_h5ad,
    ".parquet": load_parquet,
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

SUPPORTED_SUFFIXES = [sfx for sfx in FILE_LOADERS.keys() if sfx != ".rds"]
"""Suffixes with defined artifact loaders."""


def load_to_memory(
    filepath: AnyPathStr, **kwargs
) -> DataFrame | ScverseDataStructures | dict[str, Any] | list[Any] | AnyPathStr | None:
    """Load a file into memory.

    Returns the filepath if no in-memory form is found.
    May return None in interactive sessions for images.
    """
    filepath = create_path(filepath)
    suffix = extract_suffix_from_path(filepath)
    loader = FILE_LOADERS.get(suffix, None)
    if loader is None:
        raise NotImplementedError(
            f"There is no loader for {suffix} files. Use .cache() to get the path."
        )

    filepath = setup_settings.paths.cloud_to_local(filepath, print_progress=True)

    return cast(Callable[..., Any], loader)(filepath, **kwargs)
