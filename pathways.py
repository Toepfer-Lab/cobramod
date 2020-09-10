#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import logging
from typing import Union, Generator, Iterable
from creation import get_xml_from_biocyc, build_reaction_from_xml,\
    check_mass_balance, DebugLog
from cobra import Model, Reaction
from itertools import chain
from contextlib import suppress

# Creating corresponding Logs
# Format
# DebugFormatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# # Handler
# DebugHandler = logging.FileHandler("debug.log", mode="a+")
# DebugHandler.setFormatter(DebugFormatter)
# # Log
# DebugLog = logging.getLogger("DebugLog")
DebugLog.setLevel(logging.DEBUG)
# TODO!!: change name for proper python conventions


def _fix_single_reaction_dict(root: ET.Element, root_dict: dict) -> dict:
    # TODO: verify
    # These are single exceptions
    if len(root_dict) == 0:
        # single_rxn = root.find(
        #     "Pathway/reaction-list/Reaction").attrib["frameid"]
        # root_dict = {single_rxn: single_rxn}
        raise NotImplementedError(
            'Path has only a reaction. Add separately')
    return root_dict


def _create_vertex_dict(root: Union[ET.Element, str], **kwargs) -> tuple:
    """Returns a dictionary from the root of a XML file, where the keys
    represent a Reaction and the value its predecessor in a pathway. i.e, the
    direction if from right to left

    :param root: root of given pathway
    :type root: Union[ET.Element, str]
    :raises TypeError: if root is not valid
    :return: dictionary with reactions and their predecessors
    :rtype: dict
    """
    root_dict = dict()
    if isinstance(root, str):
        root = get_xml_from_biocyc(identifier=root, **kwargs)
    if not isinstance(root, ET.Element):
        msg = 'root is not a valid xml file.'
        DebugLog.error(msg)
        raise TypeError(msg)
    root_set = set()
    for rxn_line in root.findall("*/reaction-ordering"):
        current = rxn_line.find("*/[@frameid]").attrib["frameid"]
        prior = rxn_line.find("predecessor-reactions/").attrib["frameid"]
        # If the direction of keys and values changes, then
        # many reactions would get lost. This way, even with
        # multiple compounds, pathways remain
        # NOTE: check for behaviour
        # Replacing values produces cuts with are needed to avoid cyclic
        # No reactions are missed
        root_dict[current] = prior
        # TODO: add information
        root_set.add(current)
        root_set.add(prior)
    name = root.find("Pathway").attrib["frameid"]
    # NOTE: Add correspoding test
    root_dict = _fix_single_reaction_dict(root=root, root_dict=root_dict)
    DebugLog.debug(f'Dictionary for pathway "{name}" succesfully created')
    return root_dict, root_set


def _find_end_vertex(vertex_dict: dict) -> Generator:
    """Yields for a dictionary of edges, the last vertex. i.e, a single value
    that is not in a key. They can be seen as end-metabolites

    :param vertex_dict: dictionary with edges
    :type vertex_dict: dict
    :yield: last vertex
    :rtype: Generator
    """
    for value in vertex_dict.values():
        if value not in vertex_dict.keys():
            yield value


def _find_start_vertex(vertex_dict: dict) -> Generator:
    """Yields the start vertex for given dictionary, i.e. keys are do not
    appear as values. They can be seen as start-metabolites

    :param vertex_dict: dictionary with edges
    :type vertex_dict: dict
    :yield: start vertex
    :rtype: Generator
    """
    for key in vertex_dict.keys():
        if key not in vertex_dict.values():
            yield key


def _verify_return_end_vertex(vertex_dict: dict) -> Generator:
    """Verifies that at least one end-vertex and start-vertex are possible
    to obtained from dictionary of edges. If no start-vertex is found, the
    graph is cyclical and a edge will be cut

    :param vertex_dict: dictionary with edges
    :type vertex_dict: dict
    :return: end vertex
    :rtype: Generator
    """
    end_vertex = list(_find_end_vertex(vertex_dict=vertex_dict))
    while len(end_vertex) == 0:
        vertex_dict.popitem()
        end_vertex = list(_find_end_vertex(vertex_dict=vertex_dict))
    return _find_end_vertex(vertex_dict=vertex_dict)


def _build_graph(start_vertex: str, vertex_dict: dict) -> Generator:
    """Build sequence for given edge-dictionary, that starts with
    given start vertex. It will stop until no keys are available or a cyclical
    graph is found

    :param start_vertex: name of edge to start with
    :type start_vertex: str
    :param vertex_dict: dictionary with edges
    :type vertex_dict: dict
    :yield: Sequence that starts with the start vertex
    :rtype: Generator
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
    """Creates graph for given vertex dictionary. For each start vertex found,
    a sequence will be created. Edges can be cut if a cyclic graph is found

    :param vertex_dict: dictionary with edges
    :type vertex_dict: dict
    :return: list with sequences
    :rtype: Iterable
    """
    # This step is necesary to reduce the lenght of dictionary
    _ = _verify_return_end_vertex(vertex_dict=vertex_dict)
    start_vertex = _find_start_vertex(vertex_dict=vertex_dict)
    # a generator with list, otherwise items get lost.
    # True order is opposite
    return (list(
        _build_graph(
            start_vertex=start, vertex_dict=vertex_dict)
            )[::-1] for start in start_vertex)


def _return_missing_edges(
        complete_edges: Iterable, graph: Iterable) -> Generator:
    """Yields missing vertex in graph.

    :param complete_edges: complete edges of graph
    :type complete_edges: Iterable
    :param graph: graph to be tested
    :type graph: Iterable
    :yield: sequence of missing vertex
    :rtype: Generator
    """
    for edge in complete_edges:
        if edge not in graph:
            yield edge


def _replace_item(
        iterable: Iterable, replacement_dict: dict = {},
        **kwargs) -> Generator:
    """Replaces in a Iterable the item for its corresponding key in give
    dictionary. Returns a generator

    :param iterable: Sequence to modify
    :type iterable: Iterable
    :param replacement_dict: Dictionary where keys are to be replaced by
    given values, defaults to {}
    :type replacement_dict: dict, optional
    :yield: Either the original value or the new one
    :rtype: Generator
    """
    for item in iterable:
        if item in set(chain.from_iterable(replacement_dict.keys())):
            DebugLog.warning(
                f'Replacing "{item}" with "{replacement_dict[item]}"')
            yield replacement_dict[item]
        else:
            yield item


def _remove_item(
        iterable: Iterable, avoid_list: Iterable = [], **kwargs) -> Generator:
    """
    Returns Generator of items that are not in the avoid list
    """
    # TODO: add docstring
    for item in iterable:
        if item in avoid_list:
            DebugLog.warning(
                f'Avoiding root for "{item}"')
        else:
            yield item


def _return_verified_graph(
        vertex_dict: dict, vertex_set: set, **kwargs) -> list:
    """Returns list with ordered sequences for a edge dictionary. If missing
    vertex are missing in the sequences, they will be appended in the list

    :param vertex_dict: dictionary of edges
    :type vertex_dict: dict
    :return: verified list with sequences
    :rtype: list
    """
    # TODO: check for single reaction
    graph = list(get_graph(vertex_dict=vertex_dict))
    missing_edges = list(
        _return_missing_edges(
            graph=set(chain.from_iterable(graph)),
            complete_edges=vertex_set))
    graph.insert(0, missing_edges)
    # Fixing sequences from graph
    graph = [list(_replace_item(
        iterable=sequence, **kwargs)) for sequence in graph]
    graph = [list(_remove_item(
        iterable=sequence, **kwargs)) for sequence in graph]
    # Filtering NONES from .append()
    return list(filter(None, graph))


def return_graph_from_root(root: Union[str, ET.Element], **kwargs) -> list:
    """
    Returns graph from given XML root. All sequences are in order.

    Args:
        root (Union[ET.Element, str]): root of XML file or identifier for
            specific database

    Raises:
        AttributeError: If given root is invalid

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.
    Returns:
        list: sequence with reaction identifiers in order
    """

    if isinstance(root, str):
        root = get_xml_from_biocyc(identifier=root, **kwargs)
    if not isinstance(root, ET.Element):
        raise AttributeError('Given XML root is invalid.')
    vertex_dict, vertex_set = _create_vertex_dict(root=root, **kwargs)
    return _return_verified_graph(
        vertex_dict=vertex_dict, vertex_set=vertex_set, **kwargs)


def _create_reactions_for_iter(
        sequence: Iterable, **kwargs) -> Generator:
    """
    For each identifier in the sequence, a Reaction will be created. It returns
    a generator

    Args:
        sequence (Iterable): sequence with reaction identifiers for given
        database

    Keyword Arguments:
        root (Union[ET.Element, str]): root of XML file or identifier for
            specific database
        comparment (str): location of the reactions to take place. Defaults to
            cytosol "c"
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.

    Yields:
        Generator: Sequence with Reaction objects.
    """

    DebugLog.debug(f'Obtaining root for {sequence}')
    # From given list (which includes None values), retrieve only Reaction
    # Objects.
    return (build_reaction_from_xml(
        root=reaction, **kwargs) for reaction in sequence)


def _find_next_demand(
        model: Model, reaction_id: str, ignore_list: list = [],
        **kwargs) -> str:
    """
    Returns first metabolite found either in the product or reactant side
    of given reaction.
    
    Reversibility of the reaction is taken into consideration. A list with
    metabolites identifiers can be passed to ignored them.

    Args:
        model (Model): model to test
        reaction_id (str): reaction identifier for model
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore. Defaults to []

    Raises:
        TypeError: if model is invalid
        Warning: if no metabolite was found

    Returns:
        str: metabolite identifier to create a demand reaction
    """
    if not isinstance(model, Model):
        raise TypeError('Model is invalid')
    tmp_rxn = model.reactions.get_by_id(reaction_id)
    # Checking reversibility
    # left --> right
    if model.reactions.get_by_id(reaction_id).upper_bound > abs(
            model.reactions.get_by_id(reaction_id).lower_bound):
        tmp_list = [
            rxn.id for rxn in tmp_rxn.products if rxn.id not in ignore_list]
    # left <-- right
    elif model.reactions.get_by_id(reaction_id).upper_bound < abs(
            model.reactions.get_by_id(reaction_id).lower_bound):
        tmp_list = [
            rxn.id for rxn in tmp_rxn.reactants if rxn.id not in ignore_list]
    # left <--> right
    elif model.reactions.get_by_id(reaction_id).upper_bound == abs(
            model.reactions.get_by_id(reaction_id).lower_bound):
        # TODO: decide what to do (reversibility)
        # FIXME: isomers sometimes shows double demand
        tmp_list = [
            rxn.id for rxn in tmp_rxn.products if rxn.id not in ignore_list]
    if len(tmp_list) == 0:
        # Nothing found
        raise Warning('No metabolite found to become a demand')
    else:
        DebugLog.debug(
            f'Next demand selected for "{reaction_id}": "{tmp_list[0]}"')
        return tmp_list[0]


def _has_demand(
        model: Model, metabolite: str) -> bool:
    """
    Returns True if model has a demand reaction for given metabolite identifier
    """
    return any(
        ["DM_" in reaction.id for reaction in model.metabolites.get_by_id(
            metabolite).reactions])


def _remove_boundary_if_not_model(
        model: Model, metabolite: str, boundary: str):
    """
    Removes given type of boundary reaction for a specified metabolite if found
    in the model

    Args:
        model (Model): model to test
        metabolite (str): metabolite identifier in model
        boundary (str): type of boundary. Options: "sink", "demand"
    """
    type_prefix = {"sink": "SK_", "demand": "DM_"}
    if f'{type_prefix[boundary]}{metabolite}' in (
        reaction.id for reaction in model.metabolites.get_by_id(
            metabolite).reactions):
        model.remove_reactions([
            f'{type_prefix[boundary]}{metabolite}'])
        DebugLog.warning(
            f'{boundary.capitalize()} reaction for "{metabolite}" removed')



def _less_equal_than_x_rxns(
        model: Model, metabolite: str, x: int = 2) -> bool:
    """
    Returns True if given metabolite participates in less or equal than
    X reactions for given model.
    """
    reactions = model.metabolites.get_by_id(metabolite).reactions
    return len(reactions) <= x



def _fix_meta_for_boundaries(
        model: Model, metabolite: str, ignore_list: Iterable = [], **kwargs):
    """
    Verifies that given metabolite has enough sink and demand reaction. It will
    create and remove them if necessary.

    If a metabolite has a demand reaction, it will create a sink reaction if
    the total amount is equal or less than two. If no demand reaction is found
    it will create a sink reaction if the total reactions is 1. Sink reactions
    will be delete if total reaction is either more than equal 2, or 3.

    Args:
        model (Model): model to test
        metabolite (str): metabolite identifier for given metabolite.
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []

    Raises:
        Warning: If a metabolite is found in given ignore list.
    """
    DebugLog.debug(f'Checking "{metabolite}" in for sinks and demands')
    if metabolite in ignore_list:
        msg = f'Metabolite "{metabolite}" ignored'
        DebugLog.warning(msg)
        raise Warning(msg)
    # Check for to add sinks
    if _has_demand(model=model, metabolite=metabolite):
        if _less_equal_than_x_rxns(model=model, metabolite=metabolite, x=2):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite),
                type="sink")
            DebugLog.warning(f'Sink reaction created for "{metabolite}"')
        else:
            _remove_boundary_if_not_model(
                model=model, metabolite=metabolite, boundary="sink")
    else:
        if _less_equal_than_x_rxns(model=model, metabolite=metabolite, x=1):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite),
                type="sink")
            DebugLog.warning(f'Sink reaction created for "{metabolite}"')
        elif not _less_equal_than_x_rxns(
                model=model, metabolite=metabolite, x=2):
            _remove_boundary_if_not_model(
                model=model, metabolite=metabolite, boundary="sink")
    # remove demand if needed
    _remove_boundary_if_not_model(
        model=model, metabolite=metabolite, boundary="demand")


def _verify_side_sinks_for_rxn(
        model: Model, rxn_id: str, side: str, ignore_list: Iterable = []):
    """
    Checks for either the product or reactant side of a reactions, if
    participant-metabolites have enough sink reactions to produce a feasible
    answer. If necesary, it will creates or remove them.

    Args:
        model (Model): model to test.
        rxn_id (str): reaction identifier for given model.
        side (str): Side to verify. Options: "right", "left"
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []

    Raises:
        TypeError: if model is invalid
        TypeError: if ignore list is not iterable
        ValueError: if side is not in the options
    """
    if not isinstance(ignore_list, Iterable):
        raise TypeError('Ignore list is not iterable')
    if not isinstance(model, Model):
        raise TypeError('Model is invalid')
    if side == "right":
        metabolites = model.reactions.get_by_id(rxn_id).products
    elif side == "left":
        metabolites = model.reactions.get_by_id(rxn_id).reactants
    else:
        raise ValueError('Only valid options are "right" and "left"')
    for meta in metabolites:
        # create or remove
        with suppress(Warning):
            _fix_meta_for_boundaries(
                model=model, metabolite=meta.id, ignore_list=ignore_list)


# TODO: check if kwargs are needed
def verify_sinks_for_rxn(
        model: Model, rxn_id: str, ignore_list: Iterable = [], **kwargs):
    """
    Verifies, creates or removes sink reactions for metabolites found in given
    reaction if certain conditions are met.

    If a metabolite has a demand reaction, it will create a sink reaction if
    the total amount is equal or less than two. If no demand reaction is found
    it will create a sink reaction if the total reactions is 1. Sink reactions
    will be delete if total reaction is either more than equal 2, or 3.

    Args:
        model (Model): model to test
        rxn_id (str): reaction identifier for given model.
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []
    """
    # reactant side
    _verify_side_sinks_for_rxn(
        model=model, rxn_id=rxn_id, ignore_list=ignore_list, side="left")
    # product side
    _verify_side_sinks_for_rxn(
        model=model, rxn_id=rxn_id, ignore_list=ignore_list, side="right")


def _test_rxn_for_solution(
        model: Model, rxn_id: str, solution_range: tuple = (0.1, 1000),
        times: int = 0, **kwargs):
    """
    Checks if optimized objective function value of given model lies between
    a determinated range. Function is recursive and checks if sink reactions
    are enough or exceded. It creates a demand reaction for reaction for the
    test and removes it, if necesary.

    Args:
        model (Model): model to test reaction
        rxn_id (str): reaction identifier in given model
        solution_range (tuple, optional): range of solution to pass the test.
            Defaults to (0.1, 1000).
        times (int, optional): Track of recursions. Defaults to 0.

    Raises:
        ValueError: if solution is infeasible after many recursions. Depends
            from the amount of metabolites in the reaction

    """
    # TODO: add errors for wrong MODEL or rxn_id
    with suppress(KeyError):
        if rxn_id in kwargs["ignore_list"]:
            DebugLog.warning(
                f'Reaction "{rxn_id}" found in ignore list. Skipped')
            return
    if times == 0:
        DebugLog.info(f'Testing reaction "{rxn_id}"')
    # finding demand for testing
    nextDemand = _find_next_demand(model=model, reaction_id=rxn_id, **kwargs)
    with suppress(ValueError):
        model.add_boundary(model.metabolites.get_by_id(nextDemand), "demand")
        model.reactions.get_by_id(f'DM_{nextDemand}').lower_bound = 0.2
        DebugLog.debug(f'Demand "DM_{nextDemand}" added')
    # Setting maximum times for recursion
    if times == len(model.reactions.get_by_id(rxn_id).metabolites):
        msg = f'Reaction "{rxn_id}" did not passed.'
        DebugLog.critical(msg)
        raise ValueError(msg)
    # answer must be reasonable and lie between given ranges
    # comparison must be using absolute values
    if not solution_range[0] <= abs(
            model.slim_optimize()) <= solution_range[1]:
        # Append to log
        DebugLog.debug(f'Reaction "{rxn_id}" not in range')
        verify_sinks_for_rxn(
            model=model, rxn_id=rxn_id, **kwargs)
        # Recursive with 'extra' argument
        _test_rxn_for_solution(
            model=model, rxn_id=rxn_id, solution_range=solution_range,
            times=times + 1, **kwargs)
    else:
        # if works, pass and return old objective
        DebugLog.debug(f'Reaction "{rxn_id}" showed a feasible answer.')
        _remove_boundary_if_not_model(
            model=model, metabolite=nextDemand, boundary="demand")
        verify_sinks_for_rxn(
            model=model, rxn_id=rxn_id, **kwargs)


def _add_sequence(
        model: Model, sequence: list, avoid_list: Iterable = [], **kwargs):
    """
    From a sequence of Reaction objects, add each Reaction into given model. It
    checks if new reaction does not break the metabolic system.

    Args:
        model (Model): model to add reactions to.
        sequence (list): List with Reaction objects
        avoid_list (Iterable, optional): A sequence of formatted reactions to
            avoid adding to the model. This is usefull for long pathways,
            where X reactions need to be excluded. Defaults to [].

    Raises:
        TypeError: if model is invalid
        TypeError: if reactions are not valid Reaction objects
    """
    if not isinstance(model, Model):
        raise TypeError('Model is invalid')
    if not all((isinstance(reaction, Reaction) for reaction in sequence)):
        raise TypeError('Reactions are not valid objects. Check list')
    # only if not in model
    for rxn in sequence:
        if rxn.id in avoid_list:
            DebugLog.warning(
                f'Reaction "{rxn.id}" found in avoid list. Skipping.')
            continue
        if rxn.id not in model.reactions:
            model.add_reactions([rxn])
            DebugLog.info(f'Reaction "{rxn.id}" added to model')
            _test_rxn_for_solution(model=model, rxn_id=rxn.id, **kwargs)
            check_mass_balance(model=model, rxn_id=rxn.id)
        else:
            # FIXME: avoid creating reaction
            DebugLog.warning(f'Reaction "{rxn.id}" was found in model')
    DebugLog.debug('Pathway added to Model')


def _add_graph_from_root(
        model: Model, root: Union[ET.Element, str], **kwargs):
    """
    Adds root of a pathway into given model.

    Args:
        model (Model): model to append root
        root (Union[ET.Element, str]): root of XML file or identifier for
            specific database
        **kwargs: Same as 'add_graph_to_model'
    """
    # Retrieving and creating Pathway with Reactions
    # original_objective, original_direction = model.objective, \
    #     model.objective_direction
    graph = return_graph_from_root(root=root, **kwargs)
    for pathway in graph:
        sequence = list(_create_reactions_for_iter(
            sequence=pathway, model=model, **kwargs))
        _add_sequence(
            model=model, sequence=sequence, **kwargs)
    # TODO: Fix sink for metabolites in last sequence


def _add_graph_from_sequence(
        model: Model, sequence: Iterable, **kwargs):
    """
    Adds a sequence of identifiers to given model

    Args:
        model (Model): model to append graph
        sequence (Iterable): Sequence with the idenfiers to add.
        **kwargs: Same as 'add_graph_to_model'
    """
    # TODO: add tests, and docstrings
    sequence = list(_create_reactions_for_iter(
        model=model, sequence=sequence, **kwargs))
    _add_sequence(
            model=model, sequence=sequence, **kwargs)


def add_graph_to_model(
        model: Model, graph: Union[list, str, ET.Element, set], **kwargs):
    """
    Adds given graph for a pathways or a sequence of reactions identifiers
    into given model. Data will be downloaded and structured according
    to the specified database

    Args:
        model (Model): model to append graph
        graph (Union[list, str, ET.Element, set]): The identifier for the
            pathway or a iterator with the ids of the reactions

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        comparment (str): location of the reactions to take place. Defaults to
            cytosol "c"
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions.
        avoid_list (Iterable, optional): A sequence of formatted reactions to
            avoid adding to the model. This is usefull for long pathways,
            where X reactions need to be excluded.
        show_wrong (bool): For each new reaction, if it is unbalance, it will
            show the difference for the metabolites. Defaults to True
        stop_wrong (bool): For each new reaction, if it is unbalance, it will
            stop and raise an error. Defaults to False
        solution_range (tuple, optional): range of solution to pass the test
            for each new reaction. Defaults to (0.1, 1000).
    """
    # TODO: add tests, and docstrings
    if isinstance(graph, (ET.Element, str)):
        _add_graph_from_root(model=model, root=graph, **kwargs)
    elif isinstance(graph, (list, set)):
        _add_graph_from_sequence(model=model, sequence=graph, **kwargs)
