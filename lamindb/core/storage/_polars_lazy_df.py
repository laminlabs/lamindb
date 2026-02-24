from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING

from lamindb_setup.core.upath import _ensure_sync_with_fs, get_storage_region

if TYPE_CHECKING:
    from collections.abc import Iterator

    from polars import LazyFrame as PolarsLazyFrame
    from upath import UPath

POLARS_SUFFIXES = (".parquet", ".csv", ".ndjson", ".ipc")


def _polars_options(storepath: UPath) -> dict:
    polars_options: dict = {}
    storage_options: dict[str, str | bool] = {}

    fs = storepath.fs
    fs.connect()

    endpoint_url = fs.endpoint_url
    if endpoint_url is not None:
        storage_options["aws_virtual_hosted_style_request"] = False
        storage_options["aws_endpoint_url"] = endpoint_url
        if endpoint_url.startswith("http://"):
            storage_options["aws_allow_http"] = True
    else:
        storage_options["aws_region"] = get_storage_region(storepath)

    if fs.anon:
        storage_options["aws_skip_signature"] = True
    else:
        aws_key = fs.key
        aws_secret = fs.secret
        aws_token = fs.token
        if aws_key is not None and aws_secret is not None:
            storage_options["aws_access_key_id"] = aws_key
            storage_options["aws_secret_access_key"] = aws_secret
            if aws_token is not None:
                storage_options["aws_session_token"] = aws_token
        else:
            from aiobotocore.credentials import AioRefreshableCredentials

            if isinstance(
                refreshable_credentials := fs.session._credentials,
                AioRefreshableCredentials,
            ):
                refresh_sync = _ensure_sync_with_fs(
                    refreshable_credentials._refresh, fs
                )

                def credential_provider_fn():
                    # refresh and access the credentials
                    refresh_sync()
                    aws_key = refreshable_credentials._access_key
                    aws_secret = refreshable_credentials._secret_key
                    aws_token = refreshable_credentials._token
                    expiry_time = refreshable_credentials._expiry_time
                    return {
                        "aws_access_key_id": aws_key,
                        "aws_secret_access_key": aws_secret,
                        "aws_session_token": aws_token,
                    }, int(expiry_time.timestamp()) if expiry_time is not None else None

                polars_options["credential_provider"] = credential_provider_fn

    polars_options["storage_options"] = storage_options

    return polars_options


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
    if (
        not use_fsspec
        and path0.protocol == "s3"
        and "storage_options" not in kwargs
        and "credential_provider" not in kwargs
    ):
        kwargs.update(_polars_options(path0))

    open_files = []

    try:
        for path in path_list:
            open_files.append(path.open(mode="rb") if use_fsspec else path.as_posix())

        yield scans[path_list[0].suffix](open_files, **kwargs)
    finally:
        if use_fsspec:
            for open_file in open_files:
                open_file.close()
