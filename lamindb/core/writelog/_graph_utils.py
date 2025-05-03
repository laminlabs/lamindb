from typing import Optional, TypeVar

T = TypeVar("T")


def topological_sort(graph: dict[T, set[T]]) -> Optional[list[T]]:
    """Performs a topological sort on a directed acyclic graph.

    Args:
        graph: A dictionary representing an adjacency list where graph[node] is a list
              of nodes that the 'node' points to.

    Returns:
        A list of nodes in topological order, or None if a cycle is detected
    """
    # Create a dictionary to track in-degrees for each node
    in_degree = dict.fromkeys(graph, 0)

    # Calculate in-degrees for each node
    for node in graph:
        for neighbor in graph[node]:
            in_degree[neighbor] = in_degree.get(neighbor, 0) + 1

    # Find all nodes with no incoming edges
    queue = [node for node in graph if in_degree[node] == 0]
    result = []

    # Process nodes with no incoming edges
    while queue:
        node = queue.pop(0)
        result.append(node)

        # Remove this node and update in-degrees
        for neighbor in graph[node]:
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)

    # If we couldn't visit all nodes, there's a cycle
    if len(result) != len(graph):
        return None

    return result


def find_cycle(graph: dict[T, set[T]]) -> Optional[list[T]]:
    def dfs(node: T, visited: dict[T, bool], path: list[T]):
        visited[node] = True
        path.append(node)

        for neighbor in graph[node]:
            if not visited[neighbor]:
                cycle_from_neighbor = dfs(neighbor, visited, path)
                if cycle_from_neighbor is not None:
                    return cycle_from_neighbor
            elif neighbor in path:
                return path

        path.pop()
        return None

    visited = dict.fromkeys(graph.keys(), False)

    for node in graph.keys():
        if not visited[node]:
            cycle_path: list[T] = []
            if dfs(node, visited, cycle_path):
                return cycle_path + [node]
    return None
