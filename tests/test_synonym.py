import lnschema_bionty as lb

import lamindb as ln  # noqa


def test_map_synonyms():
    lb.settings.species = "human"

    # synonym not in the database
    result = lb.Gene.map_synonyms(["ABC1", "PDCD1"])
    assert result == ["HEATR6", "PDCD1"]

    mapper = lb.Gene.map_synonyms(["ABC1", "PDCD1"], return_mapper=True)
    assert mapper == {"ABC1": "HEATR6"}

    # synonym already in the database
    lb.Gene.from_bionty(symbol="BRCA2").save()
    mapper = lb.Gene.map_synonyms(["ABC1", "FANCD1"], return_mapper=True)
    assert mapper == {"FANCD1": "BRCA2", "ABC1": "HEATR6"}
    assert lb.Gene.map_synonyms(["BRCA2"]) == ["BRCA2"]
