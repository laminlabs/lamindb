from lamindb.core.writelog._graph_utils import find_cycle, topological_sort


def test_topo_sort():
    graph = {1: {2, 3}, 2: {4, 5}, 3: {}, 4: {}, 5: {}}

    topo_sort = topological_sort(graph)
    assert topo_sort is not None

    assert topo_sort.index(2) < topo_sort.index(4)
    assert topo_sort.index(2) < topo_sort.index(5)
    assert topo_sort.index(1) < topo_sort.index(2)
    assert topo_sort.index(1) < topo_sort.index(3)


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

    assert topological_sort(graph) is None

    cycle = find_cycle(graph)
    assert cycle == [1, 2, 1]


def test_find_cycle_in_tree():
    graph = {1: {2, 3}, 2: {4, 5}, 3: {}, 4: {}, 5: {}}
    assert find_cycle(graph) is None
