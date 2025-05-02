from contextlib import contextmanager

from upath import UPath

POLARS_SUFFIXES = (".parquet", ".csv", ".ndjson", ".ipc")


@contextmanager
def _open_polars_lazy_df(paths: UPath | list[UPath], **kwargs):
    try:
        import polars as pl
    except ImportError as ie:
        raise ImportError("Please install polars: pip install polars") from ie

    if isinstance(paths, UPath):
        paths = [paths]

    open_files = [path.open(mode="rb") for path in paths]

    try:
        yield pl.scan_parquet(open_files, **kwargs)
    finally:
        for open_file in open_files:
            open_file.close()
