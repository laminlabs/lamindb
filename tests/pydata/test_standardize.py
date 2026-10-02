import numpy as np
import pandas as pd
import pytest
from lamindb.models.can_curate import (
    explode_aggregated_column_to_map,
    map_synonyms,
    not_empty_none_na,
    standardize,
    to_str,
)


@pytest.fixture(scope="module")
def genes():
    gene_symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20", "GCS"]

    records = [
        {
            "symbol": "BRCA1",
            "synonyms": "PPP1R53|RNF53|FANCS|BRCC1",
            "ensembl_gene_id": "ENSG00000012048",
        },
        {
            "symbol": "A1BG",
            "synonyms": "",
            "ensembl_gene_id": "ENSG00000121410",
        },
        {
            "symbol": "BRCA2",
            "synonyms": "FAD|FAD1|BRCC2|FANCD1|FACD|FANCD|XRCC11",
            "ensembl_gene_id": "ENSG00000139618",
        },
        {
            "symbol": "A1CF",
            "synonyms": "ACF|ACF64|APOBEC1CF|ACF65|ASP",
            "ensembl_gene_id": "ENSG00000148584",
        },
        {
            "symbol": "GCLC",
            "synonyms": "GCS|BRCA1",
            "ensembl_gene_id": "ENSG00000164889",
        },
        {
            "symbol": "UGCG",
            "synonyms": "GCS",
            "ensembl_gene_id": "ENSG00000164889",
        },
        {
            "symbol": "BRCA1-1",
            "synonyms": "BRCA1-1-synonym",
            "ensembl_gene_id": "ENSG00000012048-1",
        },
        {
            "symbol": "BRCA1-1",
            "synonyms": "BRCA1-1-synonym",
            "ensembl_gene_id": "ENSG00000012048-1-1",
        },
        {
            "symbol": "BRCA1-2",
            "synonyms": "BRCA1-1-synonym",
            "ensembl_gene_id": "ENSG00000012048-2",
        },
        {
            "symbol": "Brca1",
            "synonyms": "BRCC1",
            "ensembl_gene_id": "ENSMUSG00000017146",
        },
    ]

    df = pd.DataFrame.from_records(records)

    return gene_symbols, df


def test_map_synonyms(genes):
    gene_symbols, df = genes

    mapping = map_synonyms(df=df, identifiers=gene_symbols, field="symbol")
    assert mapping == ["A1CF", "A1BG", "BRCA2", "FANCD20", "GCLC"]

    # no synonyms
    mapping = map_synonyms(df=df, identifiers=["BRCA1", "A1BG"], field="symbol")
    assert mapping == ["BRCA1", "A1BG"]


def test_map_synonyms_field_synonym(genes):
    _, df = genes

    mapping = map_synonyms(df=df, identifiers=["BRCA1", "BRCC1"], field="symbol")
    assert mapping == ["BRCA1", "BRCA1"]


def test_map_synonyms_return_mapper(genes):
    gene_symbols, df = genes

    mapper = map_synonyms(
        df=df, identifiers=gene_symbols, field="symbol", return_mapper=True
    )

    assert mapper == {"FANCD1": "BRCA2", "GCS": "GCLC"}


def test_map_synonyms_case_sensitive(genes):
    _, df = genes

    mapping = map_synonyms(
        df=df, identifiers=["A1CF", "FANCD1", "a1CF", "fancd1"], field="symbol"
    )
    assert mapping == ["A1CF", "BRCA2", "A1CF", "BRCA2"]

    mapping = map_synonyms(
        df=df,
        identifiers=["A1CF", "FANCD1", "a1CF", "fancd1"],
        field="symbol",
        case_sensitive=True,
    )
    assert mapping == ["A1CF", "BRCA2", "a1CF", "fancd1"]


def test_map_synonyms_empty_values(genes):
    _, df = genes

    result = map_synonyms(
        df=df,
        identifiers=["", " ", None, "CD3", "FANCD1"],
        field="symbol",
        return_mapper=False,
    )
    assert result == ["", " ", None, "CD3", "BRCA2"]

    mapper = map_synonyms(
        df=df,
        identifiers=["", " ", None, "CD3", "FANCD1"],
        field="symbol",
        return_mapper=True,
    )
    assert mapper == {"FANCD1": "BRCA2"}


def test_map_synonyms_keep(genes):
    _, df = genes

    assert map_synonyms(
        df, identifiers=["GCS", "A1CF"], field="symbol", keep=False
    ) == [["GCLC", "UGCG"], "A1CF"]

    assert map_synonyms(
        df, identifiers=["GCS", "A1CF"], field="symbol", keep=False, return_mapper=True
    ) == {"GCS": ["GCLC", "UGCG"]}


def test_map_synonyms_unsupported_field(genes):
    gene_symbols, df = genes
    with pytest.raises(KeyError):
        map_synonyms(df=df, identifiers=gene_symbols, field="name", return_mapper=False)
    with pytest.raises(KeyError):
        map_synonyms(
            df=df,
            identifiers=gene_symbols,
            field="symbol",
            synonyms_field="name",
            return_mapper=False,
        )
    with pytest.raises(KeyError):
        map_synonyms(
            df=df,
            identifiers=gene_symbols,
            field="symbol",
            synonyms_field="symbol",
            return_mapper=False,
        )


def test_early_mismatch():
    cell_types = {
        "name": [
            "Plasmablast",
            "conventional dendritic cell",
            "plasmablast",
        ],
        "synonyms": [
            "",
            "cDC|dendritic reticular cell|DC1|type 1 DC",
            "CD27-positive|CD38-positive|CD20-negative B cell",
        ],
    }
    df = pd.DataFrame(cell_types)

    result = standardize(
        df=df,
        identifiers=["Plasmablast", "cDC"],
        field="name",
        return_field="name",
        case_sensitive=False,
        synonyms_field="synonyms",
    )
    assert result == ["Plasmablast", "conventional dendritic cell"]


def test_map_synonyms_empty_df():
    assert (
        map_synonyms(
            df=pd.DataFrame(), identifiers=[], field="name", return_mapper=True
        )
        == {}
    )
    assert map_synonyms(df=pd.DataFrame(), identifiers=[], field="name") == []


def test_to_str():
    assert to_str(pd.Index(["A", "a", None, np.nan])).tolist() == ["a", "a", "", ""]
    assert to_str(pd.Series(["A", "a", None, np.nan])).tolist() == ["a", "a", "", ""]
    assert to_str(
        pd.Series(["A", "a", None, np.nan]), case_sensitive=True
    ).tolist() == ["A", "a", "", ""]


def test_not_empty_none_na():
    assert not_empty_none_na(["a", None, "", np.nan]).loc[0] == "a"
    assert not_empty_none_na(pd.Index(["a", None, "", np.nan])).tolist() == ["a"]
    assert not_empty_none_na(
        pd.Series(["a", None, "", np.nan], index=["1", "2", "3", "4"])
    ).to_dict() == {"1": "a"}


def test_explode_aggregated_column_to_map(genes):
    _, df = genes
    assert explode_aggregated_column_to_map(
        df, agg_col="synonyms", target_col="symbol"
    ).to_dict() == {
        "ACF": "A1CF",
        "ACF64": "A1CF",
        "ACF65": "A1CF",
        "APOBEC1CF": "A1CF",
        "ASP": "A1CF",
        "BRCA1": "GCLC",
        "BRCA1-1-synonym": "BRCA1-1",
        "BRCC1": "BRCA1",
        "BRCC2": "BRCA2",
        "FACD": "BRCA2",
        "FAD": "BRCA2",
        "FAD1": "BRCA2",
        "FANCD": "BRCA2",
        "FANCD1": "BRCA2",
        "FANCS": "BRCA1",
        "GCS": "GCLC",
        "PPP1R53": "BRCA1",
        "RNF53": "BRCA1",
        "XRCC11": "BRCA2",
    }

    assert (
        explode_aggregated_column_to_map(
            df, agg_col="synonyms", target_col="symbol", keep="last"
        ).get("GCS")
        == "UGCG"
    )
    assert explode_aggregated_column_to_map(
        df, agg_col="synonyms", target_col="symbol", keep=False
    ).get("GCS") == ["GCLC", "UGCG"]


def test_to_str_categorical_series():
    df = pd.DataFrame([np.nan, None, "a"])
    df[0] = df[0].astype("category")

    assert to_str(df[0]).tolist() == ["", "", "a"]


def test_standardize(genes):
    gene_symbols, df = genes
    assert standardize(
        df, identifiers=gene_symbols, field="symbol", return_field="ensembl_gene_id"
    ) == [
        "ENSG00000148584",
        "ENSG00000121410",
        "ENSG00000139618",
        "FANCD20",
        "ENSG00000164889",
    ]
    assert standardize(
        df,
        identifiers=gene_symbols,
        field="symbol",
        return_field="ensembl_gene_id",
        return_mapper=True,
    ) == {
        "FANCD1": "ENSG00000139618",
        "GCS": "ENSG00000164889",
        "A1BG": "ENSG00000121410",
        "A1CF": "ENSG00000148584",
    }


def test_standardize_keep(genes):
    _, df = genes
    assert standardize(
        df,
        identifiers=["A1CF", "FANCD1", "BRCA1-1"],
        field="symbol",
        return_field="ensembl_gene_id",
        return_mapper=True,
    ) == {
        "FANCD1": "ENSG00000139618",
        "A1CF": "ENSG00000148584",
        "BRCA1-1": "ENSG00000012048-1",
    }
    assert standardize(
        df,
        identifiers=["A1CF", "FANCD1", "BRCA1-1"],
        field="symbol",
        return_field="ensembl_gene_id",
        return_mapper=True,
        keep="last",
    ) == {
        "FANCD1": "ENSG00000139618",
        "A1CF": "ENSG00000148584",
        "BRCA1-1": "ENSG00000012048-1-1",
    }
    assert standardize(
        df,
        identifiers=["A1CF", "FANCD1", "BRCA1-1-synonym"],
        field="symbol",
        return_field="ensembl_gene_id",
        keep=False,
        return_mapper=True,
    ) == {
        "FANCD1": "ENSG00000139618",
        "BRCA1-1-synonym": [
            "ENSG00000012048-1",
            "ENSG00000012048-1-1",
            "ENSG00000012048-2",
        ],
        "A1CF": "ENSG00000148584",
    }
    assert standardize(
        df,
        identifiers=["A1CF", "FANCD1", "BRCA1-1"],
        field="symbol",
        return_field="ensembl_gene_id",
        keep=False,
    ) == [
        "ENSG00000148584",
        "ENSG00000139618",
        ["ENSG00000012048-1", "ENSG00000012048-1-1"],
    ]


def test_map_synonyms_field_match_first(genes):
    _, df = genes

    mapping = map_synonyms(
        df=df,
        identifiers=["Brca1"],
        field="symbol",
        synonyms_field="synonyms",
    )
    assert mapping == ["Brca1"]


def test_map_synonyms_exact_match_priority():
    """Test that exact field match takes priority over synonyms"""
    df = pd.DataFrame(
        {"symbol": ["Abca1", "ABCA1"], "synonyms": ["", "ABC1|TGD|HDLDT1"]}
    )

    # Exact match should win even though ABC1 synonym points to ABCA1
    mapping = map_synonyms(df=df, identifiers=["Abca1"], field="symbol")
    assert mapping == ["Abca1"]

    # ABC1 should map to ABCA1 via synonym
    mapping = map_synonyms(df=df, identifiers=["ABC1"], field="symbol")
    assert mapping == ["ABCA1"]


def test_map_synonyms_case_insensitive_field_before_synonym():
    """Test that case-insensitive field match takes priority over synonyms"""
    df = pd.DataFrame({"symbol": ["BRCA1", "GCLC"], "synonyms": ["", "BRCA1|GCS"]})

    # brca1 should match BRCA1 field (case-insensitive), not GCLC via synonym
    mapping = map_synonyms(
        df=df, identifiers=["brca1"], field="symbol", case_sensitive=False
    )
    assert mapping == ["BRCA1"]


def test_map_synonyms_priority_chain():
    """Test full priority: exact > case-insensitive field > synonym"""
    df = pd.DataFrame(
        {"symbol": ["Gene1", "GENE1", "GENE2"], "synonyms": ["", "", "Gene1|G1"]}
    )

    mapping = map_synonyms(
        df=df,
        identifiers=["Gene1", "gene1", "G1"],
        field="symbol",
        case_sensitive=False,
    )
    # Gene1 -> exact match to Gene1
    # gene1 -> case-insensitive match to Gene1 (first occurrence)
    # G1 -> synonym match to GENE2
    assert mapping == ["Gene1", "Gene1", "GENE2"]
