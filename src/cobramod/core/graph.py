"""Module for graph algorithm

This module creates the needed functions to find out the reaction order from
a pathway. The vertex represents the reactions and the edges symbolize the
order of the reactions. I.e. the relationship between reactions. The main
function of this module:

build_graph: From given dictionary with Parent-reaction:children-reaction,
return the corresponding non-cyclic directed graph.

"""

from __future__ import annotations

from collections import Counter
from contextlib import suppress
from itertools import chain
from typing import Any, Dict, Union

from cobramod.error import GraphKeyError


def find_missing(graph: dict[str, Union[str, tuple[str]]]):
    """
    Checks whether the given graph is missing a key.

    Args:
        graph (dict): Dictionary representing the relationships between nodes.
            A node can have several edges, which should be represented in the
            form of values.
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


def find_cycle(
    graph: dict[str, Union[str, tuple[str]]], key: str, visited: list[str]
) -> list[str]:
    visited.append(key)

    value = graph.get(key, None)

    if isinstance(value, str):
        if value in visited:
            return visited

        return find_cycle(graph, value, visited)

    elif isinstance(value, tuple):
        from_tuple: list[list[str]] = list()

        for element in value:
            # In case of a set, all values must be tested as well
            old_visited = visited.copy()

            branch = find_cycle(graph, element, old_visited)

            if branch:
                from_tuple.append(branch)

        if from_tuple:
            return from_tuple.pop()

    return []


def return_cycles(graph: dict[str, Union[str, tuple[str]]]):
    """
    Returns a nested list of cyclic paths. These paths might repeat. If the
    graph does not contain a cycle, the list is empty.

    Args:
        graph (dict): Dictionary representing the relationships between nodes.
            A node can have several edges, which should be represented in the
            form of values.

    Returns:
        List: Nested list with cyclic paths or an empty list if the graph
            does not contain a cycle.
    """
    cycles: list[list[str]] = list()

    for node in graph.keys():
        cycle = find_cycle(graph=graph, key=node, visited=[])

        if cycle:
            cycles.append(cycle)

    return cycles


def cut_cycle(graph: dict[str, Union[str, None, tuple[str]]], key: str):
    """
    Changes value of the key to None in a given dictionary. It will raise an
    error if the value is a tuple.
    """
    if not isinstance(key, str):
        raise AttributeError(
            f'Given key "{key}" cannot be cut. Probably because it is a tuple'
        )
    graph[key] = None


def cut_parents(graph: dict):
    """
    Checks if multiple parents shared a common child. If so, the graph will
    replace the values of these parents to a None and leave one of the parents
    normal.


    Args:
        graph (dict): Dictionary representing the relationships between nodes.
            A node can have several edges, which should be represented in the
            form of values.
    """
    counter = Counter(graph.values())
    # Get the keys that are not one
    values = [key for key, value in counter.items() if value != 1]
    # Get candidates
    candidates: Dict[str, Any] = dict()
    for key, value in graph.items():
        if value in values:
            try:
                candidates[value].append(key)
            except KeyError:
                candidates[value] = [key]
    # Cut all candidates but the first item
    for values in candidates.values():
        for index, key in enumerate(values):
            if index == 0:
                continue
            graph[key] = None


def back(
    graph: dict[str, Union[str, tuple[str]]],
    value: str,
    path: list[str],
    stop_list: list = [],
) -> list[str]:
    """
    Returns a list with a linear path. The function creates a sequence from
    the graph dictionary until the given value is found in the sequence or in
    the given stop_list. The function does not work with graphs that contain
    some kind of cycle and will raise a recursion error.

    Args:
        graph (dict): Dictionary representing the relationships between nodes.
            A node can have several edges, which should be represented in the
            form of values.
        value (str): The value to be searched..
        path (list): The already-visited path.
        stop_list (list): Elements that trigger the function to stop, if found.

    Returns:
        List: A list with the path that ends with the specified value or with
            an element of the stop_list.

    Raises:
        RecursionError: If an element is visited more than once due to a cycle.
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


def verify_paths(paths: list, graph: dict) -> set[str]:
    """
    Returns the missing keys that are not present in the given paths list.

    Args:
        paths (list): Paths of the given graph.
        graph (dict): Dictionary representing the relationships between nodes.
            A node can have several edges, which should be represented in the
            form of values.

    Returns:
        set: Either the missing nodes or an empty set.
    """
    reactions = set(chain.from_iterable(paths))
    return set(graph.keys()).difference(reactions)


def get_paths(graph: dict, stop_list: list) -> list[list[str]]:
    """
    Returns a list with the longest paths in a graph. This only works with
    lineal directed paths.

    Args:
        graph (dict): Dictionary representing the relationships between nodes.
            A node can have several edges, which should be represented in the
            form of values.

    Returns:
        List: Paths from given graph

    Raises:
        RecursionError: if the graph is not lineal
    """
    end_nodes = [key for key, value in graph.items() if value is None]

    # Function back will raise RecursionError if not lineal
    path_generator = (
        back(graph=graph, value=node, path=[node], stop_list=stop_list)
        for node in end_nodes
    )
    paths: list[list[str]] = list()

    for _ in end_nodes:
        with suppress(StopIteration):
            single = next(path_generator)
            paths.append(single)

    # FIXME: this is a temporal solution. When multiple tuples are located in
    # the graph, back would not be able to find all member
    difference = verify_paths(paths=paths, graph=graph)
    if difference:
        # Must be independent of each other, otherwise, the order might change
        # the behavior of the visualization
        for item in difference:
            paths.append([item])

    return paths


def get_mapping(graph: dict, stop_list: list, new: list) -> list:
    """
    Gets the mapping for given graph. The mapping defines the longest paths
    for the graph without including the previous longest path. The function
    modifies given graph.

    Args:
        graph (dict): graph (dict): Dictionary representing the relationships
            between nodes. A node can have several edges, which should be
            represented in the form of values. This will be modified.
        stop_list (list): Elements that cause the function to stop when found.
        new (list): List with new mapping. Must be empty when called for
            the first time.

    Returns:
        List: A list with the new mapping. Longest path for each recursion and
            rest.
    """
    paths = get_paths(graph=graph, stop_list=stop_list)
    # Get longest and reduce graph
    try:
        longest = max(paths, key=len)
    except ValueError:
        return []
    new.append(longest)
    for item in longest:
        graph.pop(item)
    get_mapping(graph=graph, stop_list=longest, new=new)
    return new


def get_graph_dict(graph: dict[str, Union[str, tuple[str]]]) -> list[list[str]]:
    """
    Returns the mapping for the given graph. The mapping is defined as a list
    with a "core" path and its branches. Cyclic graphs will be cut to create
    a lineal direct graph

    .. note::
        If the first element of a branch is not in the "core", then it is
        in one of the branches.

    Returns:
        List: Mapping of the graph.

    Raises:
        GraphKeyError: If graph is missing a value.
    """
    # Check that all values are represented
    find_missing(graph=graph)
    # its value cannot be None. If single element then go directly to mapping
    with suppress(KeyError):
        # Cut parents if needed
        cut_parents(graph=graph)
        # Fix cycles if found
        cycles = return_cycles(graph=graph)
        # Check until no cycles are found.
        while cycles:
            # This is modify graph
            # TODO: check whether tuples are affected
            cut_cycle(graph=graph, key=cycles[0][0])  # type: ignore
            cycles = return_cycles(graph=graph)
    # This would modify the graph. Use copy
    mapping = get_mapping(graph=graph.copy(), stop_list=[], new=[])
    mapping.sort(key=len, reverse=True)
    return mapping


def filter_graph(graph: dict, avoid_list: list) -> dict:
    """
    Returns a new graph, where items are removed if found in given avoid list.
    """
    # TODO: docstrings
    new_graph = graph.copy()
    # Remove necessary keys in copy
    for reaction in avoid_list:
        with suppress(KeyError):
            del new_graph[reaction]
        for key, value in new_graph.items():
            # Avoid NoneType

            if not value:
                continue

            elif reaction in value and isinstance(value, tuple):
                new_graph[key] = tuple(
                    item for item in value if item != reaction
                )
                # Transform to str if only one element
                if len(new_graph[key]) == 1:
                    new_graph[key] = new_graph[key][0]

            elif reaction == value:
                new_graph[key] = None
    return new_graph


def build_lineal_graph(sequence: list[str]) -> dict[str, Union[str, None]]:
    """
    Creates a returns a simple lineal directed graph from given sequence
    """
    graph: dict[str, Union[None, str]] = dict()
    for index, reaction in enumerate(sequence):
        try:
            parent = reaction
            child = sequence[index + 1]
            graph[parent] = child

        except IndexError:
            # It must be the first one
            graph[reaction] = None

    return graph
