from contextlib import contextmanager

from upath import UPath

POLARS_SUFFIXES = (".parquet", ".csv", ".ndjson", ".ipc")


@contextmanager
def _open_polars_lazy_df(paths: UPath | list[UPath], **kwargs):
    try:
        import polars as pl
    except ImportError as ie:
        raise ImportError("Please install polars: pip install polars") from ie

    if isinstance(paths, list):
        path_list = paths
    elif paths.is_dir():
        path_list = [path for path in paths.rglob("*") if path.suffix != ""]
    else:
        path_list = [paths]

    scans = {
        ".parquet": pl.scan_parquet,
        ".csv": pl.scan_csv,
        ".ndjson": pl.scan_ndjson,
        ".ipc": pl.scan_ipc,
    }

    open_files = []

    try:
        for path in path_list:
            open_files.append(path.open(mode="rb"))

        yield scans[path_list[0].suffix](open_files, **kwargs)
    finally:
        for open_file in open_files:
            open_file.close()
