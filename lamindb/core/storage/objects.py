from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeAlias

from lamindb.core._compat import (
    is_polars_dataframe,
    with_package_obj,
)

if TYPE_CHECKING:
    from pandas import DataFrame
    from polars import DataFrame as PolarsDataFrame
    from polars import LazyFrame as PolarsLazyFrame

    from .types import ScverseDataStructures

    SupportedDataTypes: TypeAlias = (
        DataFrame | PolarsDataFrame | PolarsLazyFrame | ScverseDataStructures
    )
else:
    SupportedDataTypes: TypeAlias = Any


def infer_suffix(
    dmem: SupportedDataTypes, format: str | dict[str, Any] | None = None
) -> str:
    """Infer LaminDB storage file suffix from a data object."""
    if is_polars_dataframe(dmem):
        return _infer_dataframe_suffix(format)
    has_anndata, anndata_suffix = with_package_obj(
        dmem,
        "AnnData",
        "anndata",
        lambda obj: _infer_anndata_suffix(format),
    )
    if has_anndata:
        return anndata_suffix

    has_dataframe, dataframe_suffix = with_package_obj(
        dmem,
        "DataFrame",
        "pandas",
        lambda obj: _infer_dataframe_suffix(format),
    )
    if has_dataframe:
        return dataframe_suffix

    if with_package_obj(
        dmem,
        "MuData",
        "mudata",
        lambda obj: True,  # Just checking type, not calling any method
    )[0]:
        return ".h5mu"

    has_spatialdata, spatialdata_suffix = with_package_obj(
        dmem,
        "SpatialData",
        "spatialdata",
        lambda obj: _infer_spatialdata_suffix(format),
    )
    if has_spatialdata:
        return spatialdata_suffix
    else:
        raise NotImplementedError


def _infer_anndata_suffix(format: str | dict[str, Any] | None) -> str:
    assert not isinstance(format, dict)  # noqa: S101
    if format is not None:
        # should be `.h5ad`, `.`zarr`, or `.anndata.zarr`
        if format not in {"h5ad", "zarr", "anndata.zarr"}:
            raise ValueError(
                "Error when specifying AnnData storage format, it should be"
                f" 'h5ad', 'zarr', not '{format}'. Check 'format'"
                " or the suffix of 'key'."
            )
        return "." + format
    return ".h5ad"


def _infer_dataframe_suffix(format: str | dict[str, Any] | None) -> str:
    suffix = None
    if isinstance(format, str):
        suffix = format
    elif isinstance(format, dict):
        suffix = format.get("suffix")
    if suffix in {".csv", ".lance"}:
        return suffix
    return ".parquet"


def _infer_spatialdata_suffix(format: str | dict[str, Any] | None) -> str:
    if format is None:
        return ".zarr"
    if isinstance(format, str) and format in {"spatialdata.zarr", "zarr"}:
        return f".{format}"
    raise ValueError(
        "Error when specifying SpatialData storage format, it should be"
        f" 'zarr', 'spatialdata.zarr', not '{format}'. Check 'format'"
        " or the suffix of 'key'."
    )


# for types below note that local UPaths are subclasses of Path
# Path(UPath(...)) properly coerces local UPaths and throws an error for cloud UPaths


def write_to_disk(dmem: SupportedDataTypes, filepath: Path | str, **kwargs) -> None:
    """Writes the passed in memory data to disk to a specified path."""
    sorting_keys = kwargs.pop("_lamindb_sorting_columns", None)
    if sorting_keys:
        _write_sorted_parquet(dmem, filepath, sorting_keys, **kwargs)
        return
    if is_polars_dataframe(dmem):
        import polars as pl

        if isinstance(dmem, pl.LazyFrame):
            if Path(filepath).suffix == ".csv":
                dmem.sink_csv(filepath, **kwargs)
            elif Path(filepath).suffix == ".lance":
                from ._lance import write_lance_table

                write_lance_table(dmem.collect(), filepath, **kwargs)
            else:
                dmem.sink_parquet(filepath, **kwargs)
            return
        if Path(filepath).suffix == ".csv":
            dmem.write_csv(filepath, **kwargs)
        elif Path(filepath).suffix == ".lance":
            from ._lance import write_lance_table

            write_lance_table(dmem, filepath, **kwargs)
        else:
            dmem.write_parquet(filepath, **kwargs)
        return
    if with_package_obj(
        dmem,
        "AnnData",
        "anndata",
        lambda obj: _write_anndata(obj, filepath, **kwargs),
    )[0]:
        return

    if with_package_obj(
        dmem,
        "DataFrame",
        "pandas",
        lambda obj: _write_dataframe(obj, filepath, **kwargs),
    )[0]:
        return

    if with_package_obj(dmem, "MuData", "mudata", lambda obj: obj.write(filepath))[0]:
        return

    if with_package_obj(
        dmem,
        "SpatialData",
        "spatialdata",
        lambda obj: obj.write(filepath, overwrite=True),
    )[0]:
        return

    raise NotImplementedError


def _write_sorted_parquet(
    frame: SupportedDataTypes,
    filepath: Path | str,
    sorting_keys: list[tuple[str, bool]],
    **kwargs,
) -> None:
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    from lamindb.errors import InvalidArgument, ValidationError

    if Path(filepath).suffix != ".parquet":
        raise InvalidArgument("Sorting metadata requires a Parquet artifact.")
    if is_polars_dataframe(frame):
        import polars as pl

        options = dict(kwargs.pop("pyarrow_options", {}) or {})
        kwargs.pop("use_pyarrow", None)
        if kwargs.pop("maintain_order", True) is False:
            raise InvalidArgument("Sorting metadata requires maintain_order=True.")
        if "statistics" in kwargs:
            options["write_statistics"] = kwargs.pop("statistics")
        options.update(kwargs)
        options.setdefault("compression", "zstd")
        table = (
            frame.collect() if isinstance(frame, pl.LazyFrame) else frame
        ).to_arrow()
    else:
        options = kwargs
        engine = options.pop("engine", "auto")
        if engine not in {"auto", "pyarrow"}:
            raise InvalidArgument("Sorting metadata requires the PyArrow engine.")
        table = pa.Table.from_pandas(frame, preserve_index=options.pop("index", None))

    if "sorting_columns" in options:
        raise InvalidArgument("Cannot combine sorting_columns with a sorted schema.")
    ordering = [
        (name, "ascending" if ascending else "descending")
        for name, ascending in sorting_keys
    ]
    keys = table.select([name for name, _ in sorting_keys])
    # Verify the serialized values, not the input frame's category ordering.
    for i, field in enumerate(keys.schema):
        if pa.types.is_dictionary(field.type):
            keys = keys.set_column(
                i, field.name, pc.cast(keys.column(i), field.type.value_type)
            )
    indices = pc.sort_indices(keys, sort_keys=ordering)
    if (
        len(indices) > 1
        and not pc.all(
            pc.less(indices.slice(0, len(indices) - 1), indices.slice(1))
        ).as_py()
    ):
        raise ValidationError(
            "Cannot write Parquet sorting metadata: data is not sorted as declared."
        )
    pq.write_table(
        table,
        filepath,
        sorting_columns=pq.SortingColumn.from_ordering(
            table.schema, ordering, null_placement="at_end"
        ),
        **options,
    )


def _write_anndata(dmem: Any, filepath: Path | str, **kwargs) -> None:
    suffix = Path(filepath).suffix
    if suffix == ".h5ad":
        dmem.write_h5ad(filepath, **kwargs)
        return
    elif suffix == ".zarr":
        dmem.write_zarr(filepath, **kwargs)
        return
    else:
        raise NotImplementedError


def _write_dataframe(dmem: Any, filepath: Path | str, **kwargs) -> None:
    suffix = Path(filepath).suffix
    if suffix == ".csv":
        dmem.to_csv(filepath, **kwargs)
        return
    if suffix == ".lance":
        from ._lance import write_lance_table

        write_lance_table(dmem, filepath, **kwargs)
        return
    dmem.to_parquet(filepath, **kwargs)
