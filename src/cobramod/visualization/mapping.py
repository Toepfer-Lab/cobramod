#!/usr/bin/env python3
from typing import Dict, List

from cobramod.core.graph import build_graph, get_all_values


def child_map(mapping: list, dictionary: dict) -> dict:
    """
    Returns the relation parent: children as a dictionary. A key represent the
    parent and the key their corresponding children

    Args:
        mapping (list): Paths from given graph
        dictionary (dict): Original dictionary of directed graph.

    Returns:
        Dict: Relation of parent and children.
    """
    relation: Dict[str, list] = dict()
    for index, path in enumerate(mapping):
        if len(path) == 1:
            continue
        # Check for each path
        path_2: list
        for index_2, path_2 in enumerate(mapping):
            if path_2 == path:
                continue
            # check if the first item is in the original path to check
            values = get_all_values(dictionary=dictionary, keys=path)
            if path_2[0] in values:
                try:
                    relation[str(index)].append(str(index_2))
                except KeyError:
                    relation[str(index)] = [str(index_2)]
    return relation


def get_index(dictionary: dict, path: list, value) -> int:
    """
    Return index of item in a path, whose value is found.

    Args:
        dictionary (dict): Original dictionary of directed graph
        path (list): A list with the name of nodes to be searched.
        value (Any): Value to find in path

    Returns:
        int: index of path.

    Raises:
        Warning: If value is not found.
    """
    for index, item in enumerate(path):
        values = get_all_values(dictionary=dictionary, keys=[item])
        if value in values:
            return index
    # TODO: add error
    raise Warning(f'Value "{value}" not found.')


def unformatted_matrix(graph: dict) -> List[list]:
    """
    Returns an unformatted matrix from a graph. The matrix represent the
    locations of the nodes and their relationships. Graph will be cut if
    cyclic. The matrix have 0 in empty positions.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.
    Returns
        List[list]: Representation of matrix
    """
    # If graph cyclic, it will be modified
    mapping = build_graph(graph=graph)
    # By default firsts item
    longest = mapping[0]
    # TODO: change defaults of height
    relation = child_map(mapping=mapping, dictionary=graph)
    matrix = [[0] * len(longest)]
    matrix[0] = longest
    # dictionary for start positions
    start_position = {"0": 0}
    key: list
    for index, keys in relation.items():
        # child = item, parent = index
        item: str
        for item in keys:
            # Should not raised an error since we have the relation
            position_j = get_index(
                dictionary=graph,
                path=mapping[int(index)],
                # Always look for the first item
                value=mapping[int(item)][0],
            )
            # Use the start position from parent or define new start
            try:
                start = start_position[str(index)]
            except KeyError:
                start = 0
            # add 0's and extend
            start_position[item] = start + position_j + 1
            row = [0] * start_position[item]
            row.extend(mapping[int(item)])
            matrix.append(row)
    return matrix


def format_matrix(matrix: List[list]) -> List[list]:
    return matrix
