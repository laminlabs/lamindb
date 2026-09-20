"""Storage conventions for `SpatialData`, declared in ``schema.formats.zarr``.

Each key is optional and enables one check. Keys re-use the vocabulary of the
specification they check: ``zarr_format`` and ``chunk_shape`` from the zarr v3
metadata, ``multiscales`` from OME-NGFF, ``layers`` from the anndata
``encoding-type``::

    schema.formats.zarr = {
        "zarr_format": 3,
        "chunk_shape": {"y": [256, 512], "x": [256, 512]},
        "multiscales": {"scale": 2, "axes": ["y", "x"]},
        "layers": {"csc": "csc_matrix"},
    }
    schema.save()
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..errors import ValidationError

if TYPE_CHECKING:
    try:
        from spatialdata import SpatialData
    except Exception:  # pragma: no cover
        SpatialData = Any

    from .core import SpatialDataCurator
else:
    SpatialData = Any

# downscaling floors the extent, so consecutive scales never divide exactly
SCALE_RTOL = 0.05


def _scales(image: Any) -> list[tuple[str, Any]]:
    """Return ``(scale_name, DataArray)`` pairs ordered from finest to coarsest."""
    from xarray import DataTree

    if not isinstance(image, DataTree):
        return [("", image)]

    def order(name: str) -> tuple[int, str]:
        digits = "".join(filter(str.isdigit, name))
        return (int(digits) if digits else 0, name)

    return [
        (name, next(iter(image[name].ds.data_vars.values())))
        for name in sorted(image.children, key=order)
    ]


def _chunks(array: Any) -> dict[str, int]:
    """Chunk extent per dimension, empty for a non-dask array."""
    return dict(zip(array.dims, getattr(array.data, "chunksize", ())))


def _check_zarr_format(curator: SpatialDataCurator, required: int) -> list[str]:
    import zarr

    # an in-memory object has no store to inspect
    path = getattr(curator._dataset, "path", None) or getattr(
        curator._artifact, "path", None
    )
    if path is None:
        return []
    found = zarr.open_group(str(path), mode="r").metadata.zarr_format
    if found != required:
        return [
            f"Store at '{path}' is zarr v{found}, expected v{required}"
            f"\n    → re-write the object on a zarr>={required} installation"
        ]
    return []


def _check_multiscales(sdata: SpatialData, config: dict[str, Any]) -> list[str]:
    import numpy as np

    scale = float(config.get("scale", 2))
    axes = config.get("axes", ("y", "x"))
    errors = []
    for key, image in sdata.images.items():
        scales = _scales(image)
        if len(scales) < 2:
            errors.append(
                f"Image '{key}' is single-scale"
                "\n    → build a pyramid via: "
                f"Image2DModel.parse(..., scale_factors=[{scale:g}, {scale:g}])"
            )
            continue
        for (finer_name, finer), (coarser_name, coarser) in zip(scales, scales[1:]):
            for axis in axes:
                if axis not in finer.sizes or axis not in coarser.sizes:
                    continue
                ratio = finer.sizes[axis] / coarser.sizes[axis]
                if not np.isclose(ratio, scale, rtol=SCALE_RTOL, atol=0):
                    errors.append(
                        f"Image '{key}': '{finer_name}' to '{coarser_name}' downscales "
                        f"'{axis}' by {ratio:.2f}x, expected {scale:g}x"
                    )
        if not config.get("coarsest_single_chunk", True):
            continue
        name, coarsest = scales[-1]
        chunks = _chunks(coarsest)
        for axis in axes:
            if coarsest.sizes.get(axis, 0) > chunks.get(axis, float("inf")):
                errors.append(
                    f"Image '{key}': coarsest scale '{name}' spans "
                    f"{coarsest.sizes[axis]} along '{axis}' in chunks of "
                    f"{chunks[axis]}"
                    "\n    → add scales so the coarsest one is a single chunk"
                )
    return errors


def _check_chunk_shape(sdata: SpatialData, chunk_shape: dict[str, Any]) -> list[str]:
    errors = []
    for key, image in sdata.images.items():
        for scale_name, array in _scales(image):
            chunks = _chunks(array)
            for axis, bounds in chunk_shape.items():
                if axis not in chunks:
                    continue
                smallest, largest = (
                    (bounds, bounds) if isinstance(bounds, int) else bounds
                )
                chunk, extent = chunks[axis], array.sizes[axis]
                # an axis smaller than the smallest tile may be a single chunk
                if chunk == extent and extent < smallest:
                    continue
                if not smallest <= chunk <= largest:
                    label = f"{key}/{scale_name}" if scale_name else key
                    errors.append(
                        f"Image '{label}': chunk of {chunk} along '{axis}' (extent "
                        f"{extent}), expected between {smallest} and {largest}"
                        "\n    → re-chunk on parse via: "
                        f"Image2DModel.parse(..., chunks=(1, {smallest}, {smallest}))"
                    )
    return errors


def _check_layers(sdata: SpatialData, required: dict[str, str]) -> list[str]:
    from scipy.sparse import issparse

    errors = []
    for key, table in sdata.tables.items():
        for layer, encoding_type in required.items():
            if layer not in table.layers:
                errors.append(
                    f"Table '{key}' has no layer '{layer}'"
                    f"\n    → add it via: table.layers['{layer}'] = "
                    "scipy.sparse.csc_matrix(table.X)"
                )
                continue
            matrix = table.layers[layer]
            found = f"{matrix.format}_matrix" if issparse(matrix) else "array"
            if found != encoding_type:
                errors.append(
                    f"Table '{key}': layer '{layer}' is '{found}', "
                    f"expected '{encoding_type}'"
                )
    return errors


def validate_zarr_conventions(
    curator: SpatialDataCurator, spec: dict[str, Any]
) -> None:
    """Check the curator's `SpatialData` against the conventions enabled in ``spec``."""
    sdata = curator._dataset
    errors: list[str] = []
    if "zarr_format" in spec:
        errors += _check_zarr_format(curator, spec["zarr_format"])
    if "multiscales" in spec:
        errors += _check_multiscales(sdata, spec["multiscales"])
    if "chunk_shape" in spec:
        errors += _check_chunk_shape(sdata, spec["chunk_shape"])
    if "layers" in spec:
        errors += _check_layers(sdata, spec["layers"])
    if errors:
        raise ValidationError(
            "SpatialData object violates the zarr conventions of its schema:\n"
            + "\n".join(errors)
        )
