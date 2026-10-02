import numpy as np
import pandas as pd
import pytest
from lamindb.models.can_curate import inspect, validate


@pytest.fixture(scope="module")
def genes():
    data = {
        "gene symbol": ["A1CF", "A1BG", "FANCD1", "corrupted"],
        "hgnc id": ["HGNC:24086", "HGNC:5", "HGNC:1101", "corrupted"],
        "ensembl_gene_id": [
            "ENSG00000148584",
            "ENSG00000121410",
            "ENSG00000188389",
            "ENSG0000corrupted",
        ],
    }
    data = pd.DataFrame(data).set_index("ensembl_gene_id")

    records = [
        {
            "symbol": "A1BG",
            "hgnc_id": "HGNC:5",
            "ensembl_gene_id": "ENSG00000121410",
            "synonyms": "",
        },
        {
            "symbol": "BRCA2",
            "hgnc_id": "HGNC:1101",
            "ensembl_gene_id": "ENSG00000188389",
            "synonyms": "FAD|FAD1|BRCC2|FANCD1|FACD|FANCD|XRCC11",
        },
        {
            "symbol": "A1CF",
            "hgnc_id": "HGNC:24086",
            "ensembl_gene_id": "ENSG00000148584",
            "synonyms": "ACF|ACF64|APOBEC1CF|ACF65|ASP",
        },
    ]
    df = pd.DataFrame.from_records(records)

    return df, data


def test_inspect_iterable(genes):
    df, data = genes

    result = inspect(df=df, identifiers=data.index, field="ensembl_gene_id")
    assert result.validated == [
        "ENSG00000148584",
        "ENSG00000121410",
        "ENSG00000188389",
    ]
    assert result.non_validated == ["ENSG0000corrupted"]
    assert result["validated"] == [
        "ENSG00000148584",
        "ENSG00000121410",
        "ENSG00000188389",
    ]
    assert result["non_validated"] == ["ENSG0000corrupted"]
    assert result["mapped"] == [
        "ENSG00000148584",
        "ENSG00000121410",
        "ENSG00000188389",
    ]
    assert result["not_mapped"] == ["ENSG0000corrupted"]
    with pytest.raises(KeyError):
        result["unmapped"]

    result = inspect(df=df, identifiers=data["hgnc id"], field="hgnc_id")
    assert result["validated"] == ["HGNC:24086", "HGNC:5", "HGNC:1101"]
    assert result["non_validated"] == ["corrupted"]


def test_inspect_inspect_synonyms(genes):
    df, data = genes

    result = inspect(df=df, identifiers=data["gene symbol"], field="symbol")
    assert result.validated == ["A1CF", "A1BG"]
    assert result.non_validated == ["FANCD1", "corrupted"]

    result = inspect(
        df=df, identifiers=data["gene symbol"], field="symbol", inspect_synonyms=False
    )
    assert result.validated == ["A1CF", "A1BG"]
    assert result.non_validated == ["FANCD1", "corrupted"]

    df = df.drop(columns=["synonyms"])
    result = inspect(df=df, identifiers=data["gene symbol"], field="symbol")


def test_inspect_return_df(genes):
    df, data = genes

    result = inspect(
        df=df, identifiers=data.index, field="ensembl_gene_id", return_df=True
    )

    expected_df = pd.DataFrame(
        index=[
            "ENSG00000148584",
            "ENSG00000121410",
            "ENSG00000188389",
            "ENSG0000corrupted",
        ],
        data={
            "__validated__": [True, True, True, False],
        },
    )

    assert result.equals(expected_df)

    result = inspect(df=df, identifiers=data.index, field="ensembl_gene_id")
    assert result.df.equals(expected_df)


def test_inspect_empty_dup_input(genes):
    df, _ = genes

    result = inspect(
        df=df,
        identifiers=pd.Series(["A1CF", "A1BG", "A1BG", "", None, np.nan]),
        field="symbol",
    )
    assert result.validated == ["A1CF", "A1BG"]
    assert result.non_validated == []


def test_inspect_zero_identifiers():
    result = inspect(
        df=pd.DataFrame(),
        identifiers=pd.Series([]),
        field="symbol",
    )
    assert result.validated == []
    assert result.non_validated == []


def test_inspect_empty_df():
    import pandas as pd

    result = inspect(
        df=pd.DataFrame(),
        identifiers=pd.Series(["A1CF", "A1BG", "A1BG", "", None, np.nan]),
        field="symbol",
    )

    assert result.validated == []
    assert result.non_validated == ["A1CF", "A1BG"]
    assert result.frac_validated == 0

    result = inspect(
        df=pd.DataFrame(),
        identifiers=pd.Series(["A1CF", "A1BG", "A1BG", "", None, np.nan]),
        field="symbol",
        return_df=True,
        logging=False,
    )

    expected_df = pd.DataFrame(
        index=["A1CF", "A1BG", "A1BG", "", None, np.nan],
        data={
            "__validated__": [False, False, False, False, False, False],
        },
    )

    assert result.equals(expected_df)


def test_inspect_casing(genes):
    df, _ = genes
    result = inspect(
        df=df,
        identifiers=pd.Series(["a1cf", "A1BG"]),
        field="symbol",
    )
    assert result.validated == ["A1BG"]
    assert result.non_validated == ["a1cf"]


@pytest.mark.parametrize(
    "identifiers, field_values",
    [
        ([0, 1, 2], ["a", "b", "c"]),
        ([0.0, 0.1, 0.2], [str(val) for val in [0.0, 0.1, 0.2]]),
    ],
)
def test_type_compatibility_raises_with_df(identifiers, field_values):
    df = pd.DataFrame({"field_col": field_values})
    with pytest.raises(TypeError, match="Type mismatch"):
        _ = inspect(df=df, identifiers=pd.Series(identifiers), field="field_col")


def test_validate(genes):
    df, _ = genes
    assert validate(
        identifiers=["A1CF", "a1cf"],
        field_values=df["symbol"],
        case_sensitive=True,
    ).tolist() == [True, False]
    assert validate(
        identifiers=["A1CF", "a1cf"],
        field_values=df["symbol"],
        case_sensitive=False,
        return_df=True,
    ).tolist() == [True, True]
    assert validate(
        identifiers=df["symbol"],
        field_values=df["symbol"],
    ).tolist() == [True, True, True]
