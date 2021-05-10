#!/usr/bin/env python3
"""Module for graph algorithm

This module creates the needed functions to find out the reaction order from
a pathway. The vertex represent the reactions and the edges symbolize the order
of the reactions. i.e the relationship between reactions. The main function
of this module:

build_graph: From given dictionary with Parent-reaction:children-reaction,
return the corresponding non-cyclic directed graph.
"""
from contextlib import suppress
from collections import Counter
from itertools import chain
from typing import Optional, Dict, Any
from pathlib import Path

from cobra import Model

from cobramod.error import GraphKeyError
from cobramod.core.creation import get_data, _fix_name
from cobramod.utils import _first_item
from cobramod.error import NoIntersectFound


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
        from_tuple = list()
        for single in key:
            # In case of a set, all values must be tested as well
            old_visited = visited.copy()
            old_visited.append(single)
            branch = find_cycle(graph=graph, key=single, visited=old_visited)
            # Only append if not False
            if branch:
                if isinstance(branch[0], list):
                    from_tuple.extend(branch)
                    continue
                # In case of a list
                from_tuple.append(branch)
        # Return False for an empty list
        if not from_tuple:
            return False
        # Return first object
        elif len(from_tuple) == 1:
            return from_tuple[0]
        return from_tuple


def return_cycles(graph: dict):
    """
    Returns a nested list of cyclic paths. These paths might repeat. If path
    is lineal then an empty list is returned.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.

    Returns:
        List: Nested list with cyclic paths or empty list if lineal
    """
    cycles = list()  # type: ignore
    for node in graph.keys():
        cycle = find_cycle(graph=graph, key=node, visited=[])
        if cycle:
            if isinstance(cycle[0], list):
                cycles.extend(cycle)
                continue
            cycles.append(cycle)
    return cycles


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


def cut_parents(graph: dict):
    """
    Checks if multiple parens shared a common child. If so, the graph will
    replace the values of these parents to a None and leave one of the parent
    normal.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.
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


def verify_paths(paths: list, graph: dict) -> set:
    """
    Returns the missing keys that are not present in given paths list.

    Args:
        paths (list): Paths from given graph
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.

    Returns:
        set: Either the missing nodes or an empty set
    """
    reactions = set(chain.from_iterable(paths))
    return set(graph.keys()).difference(reactions)


def get_paths(graph: dict, stop_list: list) -> list:
    """
    Returns a list with the longest path in a graph. This only works with
    lineal directed paths.

    Args:
        graph (dict): Dictionary with relationship between nodes. A node can
            have multiple edges, which should be presented as values.

    Returns:
        List: Paths from given graph

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
    # FIXME: this is a temporal solution. When multiple tuples are located in
    # the graph, back would not be able to find all member
    difference = verify_paths(paths=paths, graph=graph)
    if difference:
        # Must be independent of each other, otherwise, the order might change
        # the behavior of the visualization
        for item in difference:
            paths.append([item])
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
            cut_cycle(graph=graph, key=cycles[0][0])
            cycles = return_cycles(graph=graph)
    # This would modify the graph. Use copy
    mapping = get_mapping(graph=graph.copy(), stop_list=[], new=[])
    mapping.sort(key=len, reverse=True)
    return mapping


def _fix_graph(graph: dict, avoid_list: list, replacement: dict) -> dict:
    """
    Returns a new graph, where items are replaced if found in given
    dictionary or removed if found in given avoid list.
    """
    new_graph = graph.copy()
    # Remove necessary keys in copy
    for key in graph.keys():
        if key in avoid_list:
            del new_graph[key]
    # use new copy and modify
    graph = new_graph.copy()
    # Use replacements
    for key in graph.keys():
        if key not in replacement.keys():
            continue
        new_key = replacement[key]
        for key2, value in graph.items():
            if not value:
                # In case of None
                continue
            # If tuple and found
            elif isinstance(value, tuple) and key in value:
                new_value = tuple(item.replace(key, new_key) for item in value)
            # If string and found
            elif key in value:
                new_value = value.replace(key, new_key)
            # Skip and/or update if necessary
            else:
                continue
            new_graph[key2] = new_value
        # Update key and remove old
        new_graph[new_key] = new_graph[key]
        del new_graph[key]
    return new_graph


def _format_graph(
    graph: dict,
    model: Model,
    compartment: str,
    directory: Path,
    database: str,
    model_id: str,
    avoid_list: list,
    replacement: dict,
    genome: Optional[str],
) -> dict:
    """
    Returns new formatted graph. If an item if found in given replacement dict,
    the value will be replaced; if found in avoid list, the graph will not
    contain that value. Function :func`cobramod.get_data` will use passed
    arguments to format the items. An KeyError will be raised if the format
    was not correct.
    """
    graph = _fix_graph(
        graph=graph, avoid_list=avoid_list, replacement=replacement
    )
    new_graph = graph.copy()
    # Format tuples, single strings and kes
    for key in graph.keys():
        # Get synonym/ Format key
        # Get data from files
        key_dict = get_data(
            directory=directory,
            database=database,
            identifier=key,
            model_id=model_id,
            genome=genome,
        )
        try:
            # Find synonym
            new_identifier = _first_item(
                first=model.reactions, second=key_dict["XREF"], revert=True
            )
            new_key = f"{_fix_name(name=new_identifier)}_{compartment}"
        except (NoIntersectFound, KeyError):
            # Regular transformation
            new_key = f'{key.replace("-", "_")}_{compartment}'
        # Change in new graph the values where key is found
        for key2, value in new_graph.items():
            if not value:
                # In case of None
                continue
            # If tuple and found
            elif isinstance(value, tuple) and key in value:
                new_value = tuple(item.replace(key, new_key) for item in value)
            # If string and found
            elif key in value:
                new_value = value.replace(key, new_key)
            # Skip and/or update if necessary
            else:
                continue
            new_graph[key2] = new_value
        # Change the key and remove old one
        try:
            new_graph[new_key] = new_graph[key]
        except KeyError:
            pass
        del new_graph[key]
    # Test all names
    for key in new_graph.keys():
        model.reactions.get_by_id(key)
    return new_graph


def _create_quick_graph(sequence: list) -> dict:
    """
    Creates a returns a simple lineal directed graph from given sequence. This
    function is used for adding pathways
    """
    graph = dict()
    for index, reaction in enumerate(sequence):
        try:
            parent = reaction
            child = sequence[index + 1]
            graph[parent] = child
        except IndexError:
            # It must be the first one
            graph[reaction] = None
    return graph
