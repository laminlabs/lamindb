from lamindb.core.writelog._graph_utils import find_cycle


def test_find_nontrivial_cycle():
    graph = {
        1: {2},
        2: {3},
        3: {1},
    }
    cycle = find_cycle(graph)

    assert cycle == [1, 2, 3, 1]


def test_find_trivial_cycle():
    graph = {1: {2}, 2: {1}}

    cycle = find_cycle(graph)
    assert cycle == [1, 2, 1]


def test_find_cycle_in_tree():
    graph = {1: {2, 3}, 2: {4, 5}, 3: {}, 4: {}, 5: {}}
    assert find_cycle(graph) is None
