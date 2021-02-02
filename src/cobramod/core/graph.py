#!/usr/bin/env python3
"""Module for graph algorithm

This module creates the needed functions to find out the reaction order from
a pathway. The vertex represent the reactions and the edges symbolize the order
of the reactions. i.e the relationship between reactions. The main function
of this module:

- return_graph_from_dict: Obtain a list from the data of a pathway.
"""
from itertools import chain
from typing import Generator, Iterable

from cobramod.utils import _replace_item, _remove_item


def _find_end_vertex(vertex_dict: dict) -> Generator:
    """
    Yields from a dictionary of edges, the last vertex. i.e, a single value
    that is not in a key. They can be seen as end-metabolites.
    """
    for value in vertex_dict.values():
        if value not in vertex_dict.keys():
            yield value


def _find_start_vertex(vertex_dict: dict) -> Generator:
    """
    Yields the start vertex for given dictionary, i.e. keys are do not
    appear as values. They can be seen as start-metabolites.
    """
    for key in vertex_dict.keys():
        if key not in vertex_dict.values():
            yield key


def _verify_return_end_vertex(vertex_dict: dict) -> Generator:
    """
    Verifies that at least one end-vertex and start-vertex are possible
    to obtained from dictionary of edges. If no start-vertex is found, the
    graph is cyclical and a edge will be cut.

    Args:
        vertex_dict (dict): dictionary with edges

    Returns:
        Generator: A generator that gives last vertex

    """
    end_vertex = list(_find_end_vertex(vertex_dict=vertex_dict))
    while len(end_vertex) == 0:
        vertex_dict.popitem()
        end_vertex = list(_find_end_vertex(vertex_dict=vertex_dict))
    return _find_end_vertex(vertex_dict=vertex_dict)


def _build_graph(start_vertex: str, vertex_dict: dict) -> Generator:
    """
    Build sequence for given edge-dictionary, that starts with
    given start vertex. It will stop until no keys are available or a cyclical
    graph is found.

    Args:
        start_vertex (str): name of the first
        vertex_dict (dict): dictionary with edges

    Yields:
        Generator: sequence that start with first edge
    """
    first_step = vertex_dict[start_vertex]
    count = 0
    while True:
        try:
            yield start_vertex
            start_vertex = vertex_dict[start_vertex]
            count += 1
        except KeyError:
            break
        if start_vertex == first_step and count > 1:
            break


def get_graph(vertex_dict: dict) -> Iterable:
    """
    Creates graph for given vertex dictionary. For each start vertex found,
    a sequence will be created. Edges can be cut if a cyclic graph is found.

    Args:
        vertex_dict (dict): dictionary with edges

    Returns:
        Iterable: List with sequences
    """
    # This step is necesary to reduce the lenght of dictionary
    _ = _verify_return_end_vertex(vertex_dict=vertex_dict)
    start_vertex = _find_start_vertex(vertex_dict=vertex_dict)
    # a generator with list, otherwise items get lost.
    # True order is opposite
    return (
        list(_build_graph(start_vertex=start, vertex_dict=vertex_dict))[::-1]
        for start in start_vertex
    )


def _return_missing_edges(
    complete_edges: Iterable, graph: Iterable
) -> Generator:
    """
    Compares a iterable of edges with a graph and yields missing vertex.
    """
    for edge in complete_edges:
        if edge not in graph:
            yield edge


def _return_verified_graph(
    vertex_dict: dict,
    vertex_set: set,
    avoid_list: Iterable = [],
    replacement: dict = {},
) -> list:
    """
    Returns list with an ordered sequence for a dictionary of edges.
    The value for each key in the dictionary should be its predecessor. If a
    vertex is missing in the sequence, they will be appended as a extra
    sequence.


    Args:
        vertex_dict (dict): dictionary with edges
        vertex_set (set): set with all members
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.
        avoid_list (Iterable, Optional): original identifiers to avoid in
            graph. Defaults to [].

    Returns:
        list: verified list of sequences
    """
    # TODO: check for single reaction
    graph = list(get_graph(vertex_dict=vertex_dict))
    missing_edges = list(
        _return_missing_edges(
            graph=set(chain.from_iterable(graph)), complete_edges=vertex_set
        )
    )
    graph.append(missing_edges)
    # Fixing sequences from graph
    graph = [
        list(_replace_item(iterable=sequence, replacement=replacement))
        for sequence in graph
    ]
    # Get rid off item in avoid list
    graph = [
        list(_remove_item(iterable=sequence, avoid_list=avoid_list))
        for sequence in graph
    ]
    # Filtering NONES from .append()
    return list(filter(None, graph))


def return_graph_from_dict(
    data_dict: dict, avoid_list: Iterable = [], replacement: dict = {}
) -> list:
    """
    Returns graph from given data. All sequences are in order.

    Args:
        data_dict (dict): pathway data saved in dictionary.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.
        avoid_list (Iterable, Optional): original identifiers to avoid in
            graph. Defaults to [].
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.

    Returns:
        list: sequence with reaction identifiers in order
    """

    return _return_verified_graph(
        vertex_dict=data_dict["PATHWAY"],
        vertex_set=data_dict["SET"],
        avoid_list=avoid_list,
        replacement=replacement,
    )


# NOTE: These new functions are meant to replace the previous algorithm.
# The first step is to built the relationship between nodes. To do this,
# a directed graph can be achieved through a dictionary. A value comes from a
# key. KEY -> VALUE.


def find_cycle(graph: dict, key, visited: list):
    """
    Returns a list with the cycle in the graph or False is graph is lineal.
    This function is recursive.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.
        key (str): Key in dictionary to start to find the cycle.
        visited (list): List with keys already visited.

    Returns:
        List: Members of the graph that are in a cycle:
        False: If the graph is lineal

    Raises:
        KeyError: If a key is missing its relationship. When this happens, it
            is probably a value missing.
    """
    try:
        # Get the value and check if is already in visited. If so, then it
        # is not lineal. Otherwise, finish until return Lineal
        value = graph[key]
        if value in visited:
            return visited
        elif value is None:
            return False
        visited.append(value)
        return find_cycle(graph=graph, key=value, visited=visited)
    except KeyError:
        # If not a tuple then, the graph is missing an edge
        if not isinstance(key, tuple):
            raise KeyError(f'Value for "{key}" is missing.')
        for single in key:
            # In case of a set, all values must be tested as well
            return find_cycle(graph=graph, key=single, visited=visited)
