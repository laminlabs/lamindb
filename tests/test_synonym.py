import lnschema_bionty as lb

import lamindb as ln  # noqa


def test_map_synonyms():
    lb.settings.species = "human"

    # synonym not in the database
    result = lb.Gene.map_synonyms(["ABC1", "PDCD1"])
    assert result == ["HEATR6", "PDCD1"]

    result = lb.Gene.map_synonyms(["ABC1", "PDCD1"], field=lb.Gene.symbol)
    assert result == ["HEATR6", "PDCD1"]

    mapper = lb.Gene.map_synonyms(["ABC1", "PDCD1"], return_mapper=True)
    assert mapper == {"ABC1": "HEATR6"}

    # synonym already in the database
    print(lb.Gene.bionty())
    print(lb.Gene.from_bionty(symbol="LMNA"))
    lb.Gene.from_bionty(symbol="LMNA").save()
    mapper = lb.Gene.map_synonyms(["ABC1", "LMN1"], return_mapper=True)
    assert mapper == {"LMN1": "LMNA", "ABC1": "HEATR6"}
    assert lb.Gene.map_synonyms(["LMNA"]) == ["LMNA"]
