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

    scans = {
        ".parquet": pl.scan_parquet,
        ".csv": pl.scan_csv,
        ".nbjson": pl.scan_nbjson,
        ".ipc": pl.scan_ipc,
    }

    open_files = []

    try:
        for path in paths:
            open_files.append(path.open(mode="rb"))

        yield scans[paths[0].suffix](open_files, **kwargs)
    finally:
        for open_file in open_files:
            open_file.close()
