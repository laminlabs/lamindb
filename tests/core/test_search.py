import bionty as bt
import lamindb as ln
import pytest


@pytest.fixture(scope="module")
def prepare_cell_type_registry():
    bt.CellType.filter().all().delete()
    records = [
        {
            "ontology_id": "CL:0000084",
            "name": "T cell",
            "synonyms": "T-cell|T-lymphocyte|T lymphocyte",
            "children": ["CL:0000798", "CL:0002420", "CL:0002419", "CL:0000789"],
        },
        {
            "ontology_id": "CL:0000236",
            "name": "B cell",
            "synonyms": "B-lymphocyte|B lymphocyte|B-cell",
            "children": ["CL:0009114", "CL:0001201"],
        },
        {
            "ontology_id": "CL:0000696",
            "name": "PP cell",
            "synonyms": "type F enteroendocrine cell",
            "children": ["CL:0002680"],
        },
        {
            "ontology_id": "CL:0002072",
            "name": "nodal myocyte",
            "synonyms": "P cell|myocytus nodalis|cardiac pacemaker cell",
            "children": ["CL:1000409", "CL:1000410"],
        },
    ]
    public_records = []
    for ref_record in records:
        record = bt.CellType.from_source(ontology_id=ref_record["ontology_id"])
        assert record.name == ref_record["name"]
        assert set(record.synonyms.split("|")) == set(ref_record["synonyms"].split("|"))
        public_records.append(record)
    ln.save(public_records)
    yield "prepared"
    bt.CellType.filter().all().delete()


def test_search_synonyms(prepare_cell_type_registry):
    result = bt.CellType.search("P cell").df()
    assert set(result.name.iloc[:2]) == {"nodal myocyte", "PP cell"}


def test_search_limit(prepare_cell_type_registry):
    result = bt.CellType.search("P cell", limit=1).df()
    assert len(result) == 1


def test_search_case_sensitive(prepare_cell_type_registry):
    result = bt.CellType.search("b cell", case_sensitive=False).df()
    assert result.name.iloc[0] == "B cell"


def test_search_None():
    with pytest.raises(
        ValueError, match="Cannot search for None value! Please pass a valid string."
    ):
        bt.CellType.search(None)
