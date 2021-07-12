#!/usr/bin/env python3
"""Mapping of graph

This modules handles the mapping of a directed graph into a matrix that can be
used in Escher. The main function is
:func:`cobramod.visualization.mapping.get_mapping`, which checks if the graph
is lineal and creates the representation in a matrix. In case of cyclic path,
it will cut it.
"""
from itertools import chain
from typing import Dict, List

from cobramod.core.graph import build_graph


def get_all_values(dictionary: dict, keys: list) -> set:
    """
    Return a set with the all the values for a given list of keys. Elements of
    tuples will be added separately.

    Args:
        dictionary (dict): dictionary with keys and values
        keys (list): List of keys to get values

    Returns:
        Set: Values from given keys
    """
    set_values = set()
    for item in keys:
        value = dictionary[item]
        # In case of tuples
        if isinstance(value, tuple):
            for item in value:
                if item is not None:
                    set_values.add(item)
        else:
            if value is not None:
                set_values.add(value)
    return set_values


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
        # TODO: check if this condition is necesary
        # if len(path) == 1:
        #     continue
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
    cyclic. Unrelated nodes will be appended separately at the end. The matrix
    have 0 in empty positions.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.
    Returns
        List[list]: Representation of matrix

    Raises:
        KeyError: If keys in graph are missing
    """
    # If graph cyclic, it will be modified
    mapping = build_graph(graph=graph)
    if not mapping:
        return [[]]
    # By default firsts item
    longest = mapping[0]
    # TODO: change defaults of height
    relation = child_map(mapping=mapping, dictionary=graph)
    # If path is completely unrelated
    if not relation:
        return mapping
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
    # Add unrelated paths
    for index, line in enumerate(mapping):
        if str(index) not in relation.keys() and str(index) not in list(
            chain.from_iterable(relation.values())
        ):
            matrix.append(line)
    return matrix


def fill_matrix(matrix: List[list], length: int):
    """
    Fills each row of the matrix with 0 until it reaches given length

    Args:
        matrix (List[list]): Matrix to fill
        length (int): Desired length.
    """
    for row in matrix:
        while len(row) < length:
            row.append(0)


def format_matrix(matrix: List[list], max_length: int) -> List[list]:
    """
    Formats given matrix and returns, if possible, a reduced matrix. Matrix
    will be filled with 0's if needed.

    Args:
        matrix (List[list]): Matrix to fill
        max_length (int): Desired length.
    """
    # TODO: check if sorting is necessary
    # Sort and fill missing 0
    matrix.sort(key=len, reverse=True)
    fill_matrix(matrix=matrix, length=max_length)
    # Squeeze rows with 0s if possible
    # new_matrix: List[list] = matrix.copy()
    for index_j, row in enumerate(matrix):
        # Find row with 0s
        if 0 not in row:
            # new_matrix.append(row)
            continue
        # Check if non-zero values can be appended above
        previous = {
            index_i
            for index_i, item, in enumerate(matrix[index_j - 1])
            if item != 0
        }
        # check if all positions above are non-zero
        positions = {index_i for index_i, item, in enumerate(row) if item != 0}
        # if position is in previous then there is no space
        if positions.issubset(previous):
            # new_matrix.append(row)
            continue
        previous_row: list = matrix[index_j - 1].copy()
        for index_i, item in enumerate(row):
            # Non-zero positions
            if index_i not in positions:
                continue
            previous_row[index_i] = row[index_i]
        # Remove previous and replace with new row
        matrix[index_j] = previous_row
        del matrix[index_j - 1]
    return matrix


def get_mapping(graph: dict) -> List[list]:
    """
    Returns a matrix for the representation of given graph.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.
    Returns
        List[list]: Representation of matrix

    Raises:
        KeyError: If keys in graph are missing
    """
    # If graph is cyclic, it will be modified. Work with copy
    matrix = unformatted_matrix(graph=graph.copy())
    # Fill with 0 and merge rows if needed
    longest = len(matrix[0])
    matrix = format_matrix(matrix=matrix, max_length=longest)
    return matrix


def transpose(matrix: List[list]) -> List[list]:
    """
    Returns a transpose matrix for given matrix.
    """
    return list(map(list, zip(*matrix)))
