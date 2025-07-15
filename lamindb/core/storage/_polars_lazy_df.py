from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING

from lamindb_setup.core.upath import get_storage_region

if TYPE_CHECKING:
    from collections.abc import Iterator

    from polars import LazyFrame as PolarsLazyFrame
    from upath import UPath

POLARS_SUFFIXES = (".parquet", ".csv", ".ndjson", ".ipc")


def _polars_storage_options(storepath: UPath) -> dict[str, str | bool]:
    storage_options: dict[str, str | bool] = {}
    s3fs_options = storepath.storage_options

    endpoint_url = s3fs_options.get("endpoint_url", None)
    if endpoint_url is not None:
        storage_options["aws_virtual_hosted_style_request"] = False
        storage_options["aws_endpoint_url"] = endpoint_url
        if endpoint_url.startswith("http://"):
            storage_options["aws_allow_http"] = True
    else:
        storage_options["aws_region"] = get_storage_region(storepath)

    if s3fs_options.get("anon", False):
        storage_options["aws_skip_signature"] = True
    else:
        if "key" in s3fs_options:
            storage_options["aws_access_key_id"] = s3fs_options["key"]
        if "secret" in s3fs_options:
            storage_options["aws_secret_access_key"] = s3fs_options["secret"]
        if "token" in s3fs_options:
            storage_options["aws_session_token"] = s3fs_options["token"]

    return storage_options


@contextmanager
def _open_polars_lazy_df(
    paths: UPath | list[UPath], use_fsspec: bool = False, **kwargs
) -> Iterator[PolarsLazyFrame]:
    try:
        import polars as pl
    except ImportError as ie:
        raise ImportError("Please install polars: pip install polars") from ie

    scans = {
        ".parquet": pl.scan_parquet,
        ".csv": pl.scan_csv,
        ".ndjson": pl.scan_ndjson,
        ".ipc": pl.scan_ipc,
    }

    path_list = []
    if isinstance(paths, Path):
        paths = [paths]
    for path in paths:
        # assume http is always a file
        if getattr(path, "protocol", None) not in {"http", "https"} and path.is_dir():
            path_list += [p for p in path.rglob("*") if p.suffix != ""]
        else:
            path_list.append(path)
    # assume the filesystem is the same for all
    # it is checked in _open_dataframe
    path0 = path_list[0]
    storage_options = None
    if not use_fsspec:
        storage_options = kwargs.pop("storage_options", None)
        if path0.protocol == "s3" and storage_options is None:
            storage_options = _polars_storage_options(path0)

    open_files = []

    try:
        for path in path_list:
            open_files.append(path.open(mode="rb") if use_fsspec else path.as_posix())

        yield scans[path_list[0].suffix](
            open_files, storage_options=storage_options, **kwargs
        )
    finally:
        if use_fsspec:
            for open_file in open_files:
                open_file.close()
