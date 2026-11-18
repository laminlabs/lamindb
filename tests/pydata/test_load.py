from pathlib import Path

import anndata as ad
import lamindb as ln
import pandas as pd
import pytest

# ruff: noqa: F811
from _dataset_fixtures import get_small_mdata, get_small_sdata  # noqa


@pytest.fixture(scope="module")
def zip_file():
    filepath = Path("test.zip")
    with open(filepath, "w") as f:
        f.write("some")
    yield filepath
    filepath.unlink()


@pytest.fixture(
    scope="module",
    params=[
        "<html><head><style>h1 {color: blue}</style></head><body><h1>Test</h1></body></html>",
        '"><script>parent.document.body.dataset.compromised = "yes"</script>'
        '<p title="a & b">Fragment</p></iframe><img src=x onerror="alert(1)">',
    ],
)
def html_filepath(request):
    filepath = Path("./tmp.html")
    with open(filepath, "w") as f:
        f.write(request.param)
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def json_filepath():
    filepath = Path("./tmp.json")
    with open(filepath, "w") as f:
        f.write('{"a": 1}')
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def csv_filepath():
    filepath = Path("./tmp.csv")
    with open(filepath, "w") as f:
        f.write("a,b\n1,2")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def tsv_filepath():
    filepath = Path("./tmp.tsv")
    with open(filepath, "w") as f:
        f.write("a\tb\n1\t2")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def parquet_filepath():
    filepath = Path("./tmp.parquet")
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    df.to_parquet(filepath)
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def yaml_filepath():
    filepath = Path("./tmp.yaml")
    with open(filepath, "w") as f:
        f.write("a: 1\nb: 2")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def image_filepath():
    filepath = Path("./tmp.png")
    with open(filepath, "w") as f:
        f.write("mock image")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def svg_filepath():
    filepath = Path("./tmp.svg")
    with open(filepath, "w") as f:
        f.write("<svg><rect width='100' height='100'/></svg>")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def rds_filepath():
    filepath = Path("./tmp.rds")
    with open(filepath, "w") as f:
        f.write("mock rds")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def local_anndata_filepath():
    return ln.examples.datasets.anndata_file_pbmc68k_test().resolve()


@pytest.fixture(scope="module")
def adata(local_anndata_filepath):
    return ad.read_h5ad(local_anndata_filepath)


def test_load_anndata(local_anndata_filepath, adata):
    artifact = ln.Artifact(local_anndata_filepath, description="test")
    assert local_anndata_filepath == artifact._local_filepath
    assert local_anndata_filepath == artifact.path
    assert local_anndata_filepath == artifact.cache()

    artifact = ln.Artifact.from_anndata(adata, description="test")
    assert artifact._memory_rep is adata
    assert artifact.load() is adata
    assert artifact._local_filepath.resolve() == artifact.cache() == artifact.path


def test_load_mudata(get_small_mdata):
    artifact = ln.Artifact.from_mudata(get_small_mdata, description="test")
    assert artifact._memory_rep is get_small_mdata
    assert artifact.load() is get_small_mdata
    assert artifact._local_filepath.resolve() == artifact.cache() == artifact.path


def test_load_spatialdata(get_small_sdata):
    artifact = ln.Artifact.from_spatialdata(get_small_sdata, description="test")
    assert artifact._memory_rep is get_small_sdata
    assert artifact.load() is get_small_sdata
    assert artifact._local_filepath.resolve() == artifact.cache() == artifact.path


def load_blobs__repr__():
    example_blobs_sdata = ln.examples.datasets.spatialdata_blobs()
    blobs_af = ln.Artifact.from_spatialdata(
        example_blobs_sdata, key="example_blobs.zarr"
    ).save()
    example_blobs_sdata = blobs_af.load()
    # Must exist and not throw errors
    assert example_blobs_sdata.__repr__


@pytest.mark.parametrize("in_notebook", [False, True])
def test_load_html(html_filepath, monkeypatch, in_notebook):
    from xml.etree import ElementTree

    import IPython.display
    from lamindb.core import loaders

    monkeypatch.setattr(loaders, "is_run_from_ipython", in_notebook)
    displayed = []
    monkeypatch.setattr(IPython.display, "display", displayed.append)
    artifact = ln.Artifact(html_filepath, key=str(html_filepath))
    result = artifact.load()
    if in_notebook:
        assert result is None
        assert len(displayed) == 1
        frame = ElementTree.fromstring(displayed[0].data)  # noqa: S314
        assert frame.tag == "iframe"
        assert frame.attrib["sandbox"] == "allow-scripts"
        assert frame.attrib["srcdoc"] == html_filepath.read_text()
        assert len(frame) == 0
    else:
        assert Path(result).resolve() == html_filepath.resolve()
        assert displayed == []


def test_load_json(json_filepath):
    artifact = ln.Artifact(json_filepath, key=str(json_filepath))
    dictionary = artifact.load()
    assert dictionary["a"] == 1


def test_no_loader(zip_file):
    artifact = ln.Artifact(zip_file, key=str(zip_file))
    with pytest.raises(NotImplementedError):
        artifact.load()


def test_load_csv(csv_filepath):
    artifact = ln.Artifact(csv_filepath, key=str(csv_filepath))
    df = artifact.load()
    assert df.iloc[0, 0] == 1
    assert df.iloc[0, 1] == 2


def test_load_tsv(tsv_filepath):
    artifact = ln.Artifact(tsv_filepath, key=str(tsv_filepath))
    df = artifact.load()
    assert df.iloc[0, 0] == 1
    assert df.iloc[0, 1] == 2


def test_load_parquet(parquet_filepath):
    artifact = ln.Artifact(parquet_filepath, key=str(parquet_filepath))
    df = artifact.load()
    assert df.iloc[0, 0] == 1
    assert df.iloc[1, 1] == 4


def test_load_yaml(yaml_filepath):
    artifact = ln.Artifact(yaml_filepath, key=str(yaml_filepath))
    data = artifact.load()
    assert data["a"] == 1
    assert data["b"] == 2


def test_load_image(image_filepath):
    artifact = ln.Artifact(image_filepath, key=str(image_filepath))
    result = artifact.load()
    assert Path(result).name == image_filepath.name


def test_load_svg(svg_filepath):
    artifact = ln.Artifact(svg_filepath, key=str(svg_filepath))
    result = artifact.load()
    assert Path(result).name == svg_filepath.name


def test_load_rds(rds_filepath, ccaplog):
    artifact = ln.Artifact(rds_filepath, key=str(rds_filepath))
    result = artifact.load()
    assert "Please use `laminr` to load `.rds` files" in ccaplog.text
    assert Path(result).name == rds_filepath.name
