from contextlib import contextmanager

from upath import UPath

POLARS_SUFFIXES = (".parquet", ".csv", ".ndjson", ".ipc")


@contextmanager
def _open_polars_lazy_df(paths: UPath | list[UPath], **kwargs):
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
    if isinstance(paths, UPath):
        paths = [paths]
    for path in paths:
        # assume http is always a file
        if path.protocol not in {"http", "https"} and path.is_dir():
            path_list += [p for p in path.rglob("*") if p.suffix != ""]
        else:
            path_list.append(path)

    open_files = []

    try:
        for path in path_list:
            open_files.append(path.open(mode="rb"))

        yield scans[path_list[0].suffix](open_files, **kwargs)
    finally:
        for open_file in open_files:
            open_file.close()
