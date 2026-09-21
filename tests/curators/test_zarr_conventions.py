import anndata as ad
import lamindb as ln
import numpy as np
import pytest
import zarr
from scipy.sparse import csc_matrix, csr_matrix
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, TableModel

RNG = np.random.default_rng(0)


@pytest.fixture
def schema() -> ln.Schema:
    # otype is not hashed, so the suffix is what keeps this schema its own record
    return ln.Schema(
        name="zarr_conventions", otype="SpatialData", itype="Feature", suffix=".zarr"
    ).save()


def image(extent: int = 1024, chunk: int = 256, scale_factors=None):
    """A `c, y, x` image, pyramidal if `scale_factors` is passed."""
    return Image2DModel.parse(
        np.zeros((1, extent, extent), dtype=np.uint8),
        dims=("c", "y", "x"),
        chunks=(1, chunk, chunk),
        scale_factors=scale_factors,
    )


def table(n_obs: int = 8, n_var: int = 4) -> ad.AnnData:
    return TableModel.parse(ad.AnnData(RNG.random((n_obs, n_var), dtype=np.float32)))


def validate(sdata: SpatialData, schema: ln.Schema, spec: dict) -> None:
    schema.formats.zarr = spec
    ln.curators.SpatialDataCurator(sdata, schema).validate()


def test_no_conventions_declared(schema: ln.Schema):
    sdata = SpatialData(images={"image": image()})
    # the hook does not fire without the key, and an empty spec enables nothing
    ln.curators.SpatialDataCurator(sdata, schema).validate()
    validate(sdata, schema, {})


def test_conventions_gate_the_validation_that_follows(schema: ln.Schema):
    curator = ln.curators.SpatialDataCurator(
        SpatialData(images={"image": image(scale_factors=[2, 2])}), schema
    )
    schema.formats.zarr = {"multiscales": {"scale": 4}}
    with pytest.raises(ln.errors.ValidationError):
        curator.validate()
    assert not curator._is_validated
    schema.formats.zarr = {"multiscales": {"scale": 2}}
    curator.validate()
    assert curator._is_validated


def test_multiscales_scale_factor(schema: ln.Schema):
    spec = {"multiscales": {"scale": 2}}
    validate(SpatialData(images={"image": image(scale_factors=[2, 2])}), schema, spec)
    # downscaling floors the extent, so the ratio is never exactly the factor
    sdata = SpatialData(images={"image": image(extent=1023, scale_factors=[2, 2])})
    validate(sdata, schema, spec)

    sdata = SpatialData(images={"image": image()})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, spec)
    assert "Image 'image' is single-scale" in error.exconly()

    sdata = SpatialData(images={"image": image(scale_factors=[4])})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, spec)
    assert "downscales 'y' by 4.00x, expected 2x" in error.exconly()
    validate(sdata, schema, {"multiscales": {"scale": 4}})


def test_multiscales_coarsest_single_chunk(schema: ln.Schema):
    # one 2x scale off 1024 leaves a coarsest scale of 512 in chunks of 256
    sdata = SpatialData(images={"image": image(scale_factors=[2])})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"multiscales": {"scale": 2}})
    assert "coarsest scale 'scale1' spans 512 along 'y'" in error.exconly()

    spec = {"multiscales": {"scale": 2, "coarsest_single_chunk": False}}
    validate(sdata, schema, spec)
    # 'z' is absent from a 2D image, so nothing is compared
    validate(sdata, schema, {"multiscales": {"scale": 2, "axes": ["z"]}})


def test_chunk_shape_bounds(schema: ln.Schema):
    sdata = SpatialData(images={"image": image(chunk=256)})
    validate(sdata, schema, {"chunk_shape": {"y": [256, 512], "x": [256, 512]}})
    validate(sdata, schema, {"chunk_shape": {"y": 256}})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"chunk_shape": {"y": 512}})
    assert "chunk of 256 along 'y' (extent 1024)" in error.exconly()

    # a single chunk smaller than the smallest tile is left alone
    small = SpatialData(images={"image": image(extent=64, chunk=64)})
    validate(small, schema, {"chunk_shape": {"y": [256, 512]}})


def test_chunk_shape_covers_every_scale(schema: ln.Schema):
    sdata = SpatialData(images={"image": image(scale_factors=[2, 2], chunk=256)})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"chunk_shape": {"y": [512, 1024]}})
    messages = error.exconly()
    assert "Image 'image/scale0'" in messages
    assert "Image 'image/scale1'" in messages
    # scale2 is a single chunk smaller than the smallest tile, so it is exempt
    assert "Image 'image/scale2'" not in messages


def test_layers_encoding(schema: ln.Schema):
    sdata = SpatialData(tables={"table": table()})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"layers": {"csc": "csc_matrix"}})
    assert "Table 'table' has no layer 'csc'" in error.exconly()

    adata = sdata.tables["table"]
    adata.layers["csc"] = csr_matrix(adata.X)
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"layers": {"csc": "csc_matrix"}})
    assert "layer 'csc' is 'csr_matrix', expected 'csc_matrix'" in error.exconly()

    adata.layers["csc"] = csc_matrix(adata.X)
    validate(sdata, schema, {"layers": {"csc": "csc_matrix"}})

    adata.layers["dense"] = np.asarray(adata.X)
    validate(sdata, schema, {"layers": {"dense": "array"}})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"layers": {"dense": "csr_matrix"}})
    assert "layer 'dense' is 'array', expected 'csr_matrix'" in error.exconly()


def test_layers_checked_per_table(schema: ln.Schema):
    with_csc, without_csc = table(), table()
    with_csc.layers["csc"] = csc_matrix(with_csc.X)
    sdata = SpatialData(tables={"with_csc": with_csc, "without_csc": without_csc})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"layers": {"csc": "csc_matrix"}})
    messages = error.exconly()
    assert "Table 'without_csc' has no layer 'csc'" in messages
    assert "Table 'with_csc'" not in messages


def test_zarr_format(schema: ln.Schema, tmp_path):
    sdata = SpatialData(images={"image": image()})
    # an in-memory object has no store to inspect
    validate(sdata, schema, {"zarr_format": 3})

    sdata.write(tmp_path / "sdata.zarr")
    found = zarr.open_group(str(sdata.path), mode="r").metadata.zarr_format
    validate(sdata, schema, {"zarr_format": found})
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, {"zarr_format": 2 if found == 3 else 3})
    assert f"is zarr v{found}" in error.exconly()


def test_complete_spec_passes(schema: ln.Schema, tmp_path):
    adata = table()
    adata.layers["csc"] = csc_matrix(adata.X)
    sdata = SpatialData(
        images={"image": image(scale_factors=[2, 2])}, tables={"table": adata}
    )
    sdata.write(tmp_path / "sdata.zarr")
    validate(
        sdata,
        schema,
        {
            "zarr_format": zarr.open_group(
                str(sdata.path), mode="r"
            ).metadata.zarr_format,
            "multiscales": {"scale": 2},
            "chunk_shape": {"y": [256, 512], "x": [256, 512]},
            "layers": {"csc": "csc_matrix"},
        },
    )


def test_violations_are_aggregated(schema: ln.Schema):
    sdata = SpatialData(images={"image": image()}, tables={"table": table()})
    spec = {
        "multiscales": {"scale": 2},
        "chunk_shape": {"y": 512},
        "layers": {"csc": "csc_matrix"},
    }
    with pytest.raises(ln.errors.ValidationError) as error:
        validate(sdata, schema, spec)
    messages = error.exconly()
    assert "Image 'image' is single-scale" in messages
    assert "chunk of 256 along 'y' (extent 1024)" in messages
    assert "Table 'table' has no layer 'csc'" in messages
