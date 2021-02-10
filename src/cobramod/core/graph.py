#!/usr/bin/env python3
"""Module for graph algorithm

This module creates the needed functions to find out the reaction order from
a pathway. The vertex represent the reactions and the edges symbolize the order
of the reactions. i.e the relationship between reactions. The main function
of this module:

- return_graph_from_dict: Obtain a list from the data of a pathway.

TODO: add new algorithm to docstrings
"""
from contextlib import suppress
from itertools import chain
from typing import Generator, Iterable

from cobramod.error import GraphKeyError
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


def find_missing(graph: dict):
    """
    Raise KeyError is graph is missing an Key

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.

    Raises:
        KeyError: If keys are missing
    """
    values = set()
    for value in graph.values():
        if isinstance(value, tuple):
            for item in value:
                values.add(item)
            continue
        values.add(value)
    keys = {key for key in graph.keys()}
    intersection = set.difference(values, keys)
    intersection = set(filter(None, intersection))
    if len(intersection) > 0:
        raise GraphKeyError(
            f"The graph is missing following keys: {intersection}"
        )


def find_cycle(graph: dict, key: str, visited: list):
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
        GraphKeyError: If a key is missing its relationship. When this happens,
            it is probably a value missing.
    """
    try:
        # Get the value and check if is already in visited. If so, then it
        # is not lineal. Otherwise, finish until return Lineal
        value = graph[key]
        if value in visited:
            return visited
        elif value is None:
            return False
        # For tuples, the value will be added in the exception
        if not isinstance(value, tuple):
            visited.append(value)
        return find_cycle(graph=graph, key=value, visited=visited)
    except KeyError:
        # If not a tuple then, the graph is missing an edge
        if not isinstance(key, tuple):
            raise GraphKeyError(f'Value for "{key}" is missing.')
        for single in key:
            # In case of a set, all values must be tested as well
            visited.append(single)
            return find_cycle(graph=graph, key=single, visited=visited)


def cut_cycle(graph: dict, key: str):
    """
    Changes value of key to None in given dictionary. It will raise an error if
    value is a tuple
    """
    if not isinstance(key, str):
        raise AttributeError(
            f'Given key "{key}" cannot be cut. Probably because it is a tuple'
        )
    graph[key] = None
    # TODO: add debug. Is it necesary?


def back(graph: dict, value: str, path: list, stop_list: list = []) -> list:
    """
    Return a list with the path that ends with given value until it reaches
    a node without a value as a key or the value can be found in the stop_list.
    The function is recursive and will only work with lineal directed graphs.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.
        value (str): The value to search.
        path (list): The already-visited path.
        stop_list (list): Elements that which trigger the function to stop,
            if found.

    Returns:
        List: a list with the path that end up with value.

    Raises:
        RecursionError: If the path is not lineal
    """
    # Check for each key if the value matches the argument value. Extend path
    # and call recursive
    for key, val in graph.items():
        if isinstance(val, tuple):
            # In case of tuples, test each single case
            for single in val:
                if single in stop_list:
                    path.insert(0, key)
                    return path
                if single == value:
                    path.insert(0, key)
                    return back(
                        graph=graph, value=key, path=path, stop_list=stop_list
                    )
        if val in stop_list:
            path.insert(0, key)
            return path
        # If val from dictionary is the same as value, then run recursive.
        if val == value:
            path.insert(0, key)
            return back(graph=graph, value=key, path=path, stop_list=stop_list)
    return path


def get_paths(graph: dict, stop_list: list) -> list:
    """
    Returns a list with the longest path in a graph. This only works with
    lineal directed paths.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.

    Returns:
        List: Longest path

    Raises:
        RecursionError: if path is not lineal
    """
    end_nodes = [key for key, value in graph.items() if value is None]
    # Function back will raise RecursionError if not lineal
    path_generator = (
        back(graph=graph, value=node, path=[node], stop_list=stop_list)
        for node in end_nodes
    )
    paths = list()
    for node in end_nodes:
        with suppress(StopIteration):
            single = next(path_generator)
            paths.append(single)
    return paths


def get_mapping(graph: dict, stop_list: list, new: list):
    """
    Gets the mapping for given graph. The mapping defines the longest paths
    for the graph without including the previous longest path. The function
    modifies given graph.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values. This will
            be modified
        stop_list (list): Elements that which trigger the function to stop,
            if found.
        new (list): List with new mapping. Should be empty when call for the
            first time.

    Returns:
        List: A list with the new mapping. Longest path for each recursion and
            rest.
    """
    paths = get_paths(graph=graph, stop_list=stop_list)
    # Get longest and reduce graph
    try:
        longest = max(paths, key=len)
    except ValueError:
        return
    new.append(longest)
    for item in longest:
        graph.pop(item)
    get_mapping(graph=graph, stop_list=longest, new=new)
    return new


def get_pop_key(dictionary: dict) -> str:
    """
    Get key from dictionary that is not None nor a tuple.
    """
    # Copy to avoid modifications
    new = dictionary.copy()
    true_value = None
    while true_value is None:
        key, value = new.popitem()
        if isinstance(value, tuple):
            for item in value:
                true_value = item
        # key, value
        else:
            true_value = value
    return key


# TODO: remove if necesary
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


def get_key(dictionary: dict, value):
    """
    Returns key for given value. Tuples analyze separately.

    Args:
        dictionary (dict): Dictionary with keys and values
        value (Any): Value to find in dictionary

    Returns:
        key: Key for specific value

    Raises:
        Warning: If value is not found
    """
    for key, test_value in dictionary.items():
        # In case of tuples
        if isinstance(test_value, tuple):
            for item in test_value:
                if item == value:
                    return key
        if test_value == value:
            return key
    raise Warning(f'Key for "{value}" not found in dictionary.')


# TODO: remove if necesary
def check_values(item: str, longest: list, mapping: list, graph: dict):
    """
    Check that given item is located in a branch. Returns the key of the item
    if found in branch that is not the longest. If nothing is found, an error
    will be raised.

    Args:
        item (str): Value to search
        longest (list): List of longest path.
        mapping (list): List with paths
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.

    Returns:
        Str: Key from value if found in a branch.

    Raises:
        Warning: If value is not found
    """
    found = False
    path: list
    for path in mapping:
        # Skip core and singles
        if longest == path or len(path) == 1:
            continue
        values = get_all_values(dictionary=graph, keys=path)
        if item in values:
            found = True
            key = get_key(dictionary=graph, value=item)
            return key
    if not found:
        # TODO: add new Error
        raise Warning(f'Value "{item}" not found in branches of mapping.')


def build_graph(graph: dict) -> list:
    """
    Returns the mapping for given graph. The mapping is the defined as a list
    with a "core" path and its branches. Cyclic graphs will be cut in order
    to create a lineal direct graph

    .. note::
        If the first item of a branch is not in the "core", the it is in one
        of the branches.

    Returns:
        List: Mapping from graph.

    Raises:
        GraphKeyError: If graph is missing a value.
    """
    # Check that all values are represented
    find_missing(graph=graph)
    # its value cannot be None. If single element then go direcly to mapping
    with suppress(KeyError):
        key = get_pop_key(dictionary=graph)
        # Fix cycles if found
        cycle = find_cycle(graph=graph, key=key, visited=[])
        if cycle is not False:
            # This is modify graph
            cut_cycle(graph=graph, key=cycle[0])
    # This would modify the graph. Use copy
    mapping = get_mapping(graph=graph.copy(), stop_list=[], new=[])

    # TODO: check if this is necesary
    # longest = max(mapping, key=len)
    # Search for rest values, which are not included in core "longest"
    # path: list
    # for path in mapping:
    #     # Skip core
    #     if longest == path or len(path) == 1:
    #         continue
    #     # Check whether in this set
    #     values = get_all_values(dictionary=graph, keys=longest)
    #     if path[0] not in values:
    #         key = check_values(
    #             item=path[0], longest=longest, mapping=mapping, graph=graph
    #         )
    #         # Get key and insert it
    #         path.insert(0, key)
    # Return a sorted list
    mapping.sort(key=len, reverse=True)
    return mapping
