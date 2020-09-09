#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import logging
from typing import Union, Generator, Iterable
from creation import get_xml_from_biocyc, build_reaction_from_xml,\
    stopAndShowMassBalance, DebugLog
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


def fix_single_reaction_dict(root: ET.Element, root_dict: dict) -> dict:
    # TODO: verify
    # These are single exceptions
    if len(root_dict) == 0:
        # single_rxn = root.find(
        #     "Pathway/reaction-list/Reaction").attrib["frameid"]
        # root_dict = {single_rxn: single_rxn}
        raise NotImplementedError(
            'Path has only a reaction. Add separately')
    return root_dict


def create_vertex_dict(root: Union[ET.Element, str], **kwargs) -> tuple:
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
        root = get_xml_from_biocyc(bioID=root, **kwargs)
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
    root_dict = fix_single_reaction_dict(root=root, root_dict=root_dict)
    DebugLog.debug(f'Dictionary for pathway "{name}" succesfully created')
    return root_dict, root_set


def find_end_vertex(vertex_dict: dict) -> Generator:
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


def find_start_vertex(vertex_dict: dict) -> Generator:
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


def verify_return_end_vertex(vertex_dict: dict) -> Generator:
    """Verifies that at least one end-vertex and start-vertex are possible
    to obtained from dictionary of edges. If no start-vertex is found, the
    graph is cyclical and a edge will be cut

    :param vertex_dict: dictionary with edges
    :type vertex_dict: dict
    :return: end vertex
    :rtype: Generator
    """
    end_vertex = list(find_end_vertex(vertex_dict=vertex_dict))
    while len(end_vertex) == 0:
        vertex_dict.popitem()
        end_vertex = list(find_end_vertex(vertex_dict=vertex_dict))
    return find_end_vertex(vertex_dict=vertex_dict)


def build_graph(start_vertex: str, vertex_dict: dict) -> Generator:
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
    _ = verify_return_end_vertex(vertex_dict=vertex_dict)
    start_vertex = find_start_vertex(vertex_dict=vertex_dict)
    # a generator with list, otherwise items get lost.
    # True order is opposite
    return (list(
        build_graph(
            start_vertex=start, vertex_dict=vertex_dict)
            )[::-1] for start in start_vertex)


def return_missing_edges(
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
    # TODO: add docstring
    for item in iterable:
        if item in avoid_list:
            DebugLog.warning(
                f'Avoiding root for "{item}"')
        else:
            yield item


def return_verified_graph(
        vertex_dict: dict, vertex_set: set, **kwargs) -> list:
    """Returns list with ordered sequences for a edge dictionary. If missing
    vertex are missing in the sequences, they will be appended in the list

    :param vertex_dict: dictionary of edges
    :type vertex_dict: dict
    :return: verified list with sequences
    :rtype: list
    """
    # Backup for comparison
    # vertex_dict_copy = vertex_dict.copy()
    # TODO: check for single reaction
    graph = list(get_graph(vertex_dict=vertex_dict))
    missing_edges = list(
        return_missing_edges(
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
    """Returns graph from given XML rool. All sequences will be included in
    order.

    :param root: name of root or XML
    :type root: Union[str, ET.Element]
    :raises AttributeError: If given root is invalid
    :return: verified list with sequences
    :rtype: list
    """
    if isinstance(root, str):
        root = get_xml_from_biocyc(bioID=root, **kwargs)
    if not isinstance(root, ET.Element):
        raise AttributeError('Given XML root is invalid.')
    vertex_dict, vertex_set = create_vertex_dict(root=root, **kwargs)
    return return_verified_graph(
        vertex_dict=vertex_dict, vertex_set=vertex_set, **kwargs)


def create_reactions_for_iter(
        sequence: Iterable, **kwargs) -> Generator:
    """For each item in a sequence, create a Reaction object
    and return. TODO: add replacement_dict!!

    :param sequence: Iterable with items, normally strings
    :type sequence: Iterable
    :param replacement_dict: !!, defaults to {}
    :type replacement_dict: dict, optional
    :yield: Sequences with Reactions objects
    :rtype: Generator
    """

    DebugLog.debug(f'Obtaining root for {sequence}')
    # From given list (which includes None values), retrieve only Reaction
    # Objects.
    return (build_reaction_from_xml(
        root=reaction, **kwargs) for reaction in sequence)


def find_next_demand(
        model: Model, check_rxn_id: str, ignore_list: list = [],
        **kwargs) -> str:
    """Returns first metabolites found either in the product or reactant side
    of given reaction. Reversibility of the reaction is taken into
    consideration. A list with metabolites IDs can be passed to be ignored

    :param model: Model where the reaction is located
    :type model: Model
    :param check_rxn_id: ID of the reaction in the model
    :type check_rxn_id: str
    :param ignore_list: list with metabolites IDs to ignore, defaults to []
    :type ignore_list: list, optional
    :raises TypeError: if Model is not valid
    :raises Warning: if no metabolite was found
    :return: ID of next metabolite to become demand
    :rtype: str
    """
    if not isinstance(model, Model):
        raise TypeError('Model is invalid')
    tmp_rxn = model.reactions.get_by_id(check_rxn_id)
    # Checking reversibility
    # left --> right
    if model.reactions.get_by_id(check_rxn_id).upper_bound > abs(
            model.reactions.get_by_id(check_rxn_id).lower_bound):
        tmp_list = [
            rxn.id for rxn in tmp_rxn.products if rxn.id not in ignore_list]
    # left <-- right
    elif model.reactions.get_by_id(check_rxn_id).upper_bound < abs(
            model.reactions.get_by_id(check_rxn_id).lower_bound):
        tmp_list = [
            rxn.id for rxn in tmp_rxn.reactants if rxn.id not in ignore_list]
    # left <--> right
    elif model.reactions.get_by_id(check_rxn_id).upper_bound == abs(
            model.reactions.get_by_id(check_rxn_id).lower_bound):
        # TODO: decide what to do (reversibility)
        # FIXME: isomers sometimes shows double demand
        tmp_list = [
            rxn.id for rxn in tmp_rxn.products if rxn.id not in ignore_list]
    if len(tmp_list) == 0:
        # Nothing found
        raise Warning('No metabolite found to become a demand')
    else:
        DebugLog.debug(
            f'Next demand selected for "{check_rxn_id}": "{tmp_list[0]}"')
        return tmp_list[0]


def has_demand(
        model: Model, metabolite: str) -> bool:
    """Returns True if model has a demand reaction for given metabolite

    :param model: model to search for reactions
    :type model: Model
    :param metabolite: metabolite ID
    :type metabolite: str
    :return: Whether demand can be found
    :rtype: bool
    """
    return any(
        ["DM_" in reaction.id for reaction in model.metabolites.get_by_id(
            metabolite).reactions])


def remove_boundary_if_not_model(
        model: Model, metabolite: str, boundary: str):
    """Removes given type of boundary reaction for a specific metabolite if
    found in model

    :param model: model to check for reactions
    :type model: Model
    :param metabolite: name of metabolite
    :type metabolite: str
    :param boundary: type of boundary, Options are 'sink' or 'demand'.
    :type boundary: str
    """
    type_prefix = {"sink": "SK_", "demand": "DM_"}
    if f'{type_prefix[boundary]}{metabolite}' in (
        reaction.id for reaction in model.metabolites.get_by_id(
            metabolite).reactions):
        model.remove_reactions([
            f'{type_prefix[boundary]}{metabolite}'])
        DebugLog.warning(
            f'{boundary.capitalize()} reaction for "{metabolite}" removed')


def more_than_two_reaction(model: Model, metabolite: str) -> bool:
    """Returns True if given metabolite participates more than in two reaction
    for given model.

    :param model: model to check for reactions
    :type model: Model
    :param metabolite: name of metabolite
    :type metabolite: str
    :return: True if more than 2 reactions.
    :rtype: bool
    """
    reactions = model.metabolites.get_by_id(metabolite).reactions
    return len(reactions) > 2


def less_equal_than_x_rxns(
        model: Model, metabolite: str, x: int = 2) -> bool:
    """Returns True if given metabolite participates in less than 2 reactions
    for given model.

    :param model: model to check for reactions
    :type model: Model
    :param metabolite: name of metabolite
    :type metabolite: str
    :return: True if less than 2 reactions
    :rtype: bool
    """
    reactions = model.metabolites.get_by_id(metabolite).reactions
    return len(reactions) <= x


def check_if_boundary(model: Model, metabolite: str) -> bool:
    """Returns True if at least one kind of boundary (with exception of
    Exchange reactions) is found for given metabolite

    :param model: model to check for reactions
    :type model: Model
    :param metabolite: name of metabolite
    :type metabolite: str
    :return: If boundary reaction is found
    :rtype: bool
    """
    list_no_exchanges = [demand.id for demand in model.demands]\
        + [sink.id for sink in model.sinks]
    return any([
        f'DM_{metabolite}' in list_no_exchanges,
        f'SK_{metabolite}' in list_no_exchanges])


def fix_meta_for_boundaries(
        model: Model, metabolite: str, ignore_list: Iterable = [], **kwargs):
    DebugLog.debug(f'Checking "{metabolite}" in for sinks and demands')
    if metabolite in ignore_list:
        msg = f'Metabolite "{metabolite}" ignored'
        DebugLog.warning(msg)
        raise Warning(msg)
    # Check for to add sinks
    if has_demand(model=model, metabolite=metabolite):
        if less_equal_than_x_rxns(model=model, metabolite=metabolite, x=2):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite),
                type="sink")
            DebugLog.warning(f'Sink reaction created for "{metabolite}"')
        else:
            remove_boundary_if_not_model(
                model=model, metabolite=metabolite, boundary="sink")
    else:
        if less_equal_than_x_rxns(model=model, metabolite=metabolite, x=1):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite),
                type="sink")
            DebugLog.warning(f'Sink reaction created for "{metabolite}"')
        elif not less_equal_than_x_rxns(
                model=model, metabolite=metabolite, x=2):
            remove_boundary_if_not_model(
                model=model, metabolite=metabolite, boundary="sink")
    remove_boundary_if_not_model(
        model=model, metabolite=metabolite, boundary="demand")


def verify_side_sinks_for_rxn(
        model: Model, rxnID: str, side: str, ignore_list: Iterable = []):
    """Checks for either the product or reactant side of a reactions, if
    participant-metabolites have enough sink reactions and, if necesary,
    creates or remove them

    :param model: model where Reactions and Metabolite are located
    :type model: Model
    :param rxnID: ID of reaction in model
    :type rxnID: str
    :param side: which side to check. Options are 'right' or 'left',
    defaults to "right"
    :type side: str, optional
    :param ignore_list: list with metabolites to be ignored, defaults to []
    :type ignore_list: list, optional
    :raises TypeError: if model is invalid
    :raises TypeError: if ignore_list is not passed as a list
    :raises ValueError: If 'side' is not properly passed
    """
    if not isinstance(ignore_list, Iterable):
        raise TypeError('Ignore list is not iterable')
    if not isinstance(model, Model):
        raise TypeError('Model is invalid')
    if side == "right":
        metabolites = model.reactions.get_by_id(rxnID).products
    elif side == "left":
        metabolites = model.reactions.get_by_id(rxnID).reactants
    else:
        raise ValueError('Only valid options are "right" and "left"')
    for meta in metabolites:
        # create or remove
        with suppress(Warning):
            fix_meta_for_boundaries(
                model=model, metabolite=meta.id, ignore_list=ignore_list)


# TODO: check if kwargs are needed
def verify_sinks_for_rxn(
        model: Model, rxnID: str, ignore_list: Iterable = [], **kwargs):
    """Checks, creates or remove sink reactions for metabolites in given
    Reaction

    :param model: model with Reactions
    :type model: Model
    :param rxnID: ID of the reaction in model
    :type rxnID: str
    :param ignore_list: metabolites to be ignored, defaults to []
    :type ignore_list: list, optional
    """
    # reactant side
    verify_side_sinks_for_rxn(
        model=model, rxnID=rxnID, ignore_list=ignore_list, side="left")
    # product side
    verify_side_sinks_for_rxn(
        model=model, rxnID=rxnID, ignore_list=ignore_list, side="right")


def test_rxn_for_solution(
        model: Model, rxnID: str, solution_range: tuple = (0.1, 1000),
        times: int = 0, **kwargs):
    """Checks if optimized objective function value of given model lies between
    a determinated range. Function is recursive and checks if sink reactions
    are enough or exceded. It creates a demand reaction for reaction for the
    test and removes it, if necesary.

    :param model: model to test
    :type model: Model
    :param rxnID: reaction ID to test
    :type rxnID: str
    :param solution_range: range of solution, defaults to (0.1, 1000)
    :type solution_range: tuple, optional
    :param times: For recursion, displays how many times the FUN called itself,
    defaults to 0
    :type times: int, optional
    :raises ValueError: if solution is infeasible after many recursions
    """
    # TODO: add errors for wrong MODEL or rxnid
    with suppress(KeyError):
        if rxnID in kwargs["ignore_list"]:
            DebugLog.warning(
                f'Reaction "{rxnID}" found in ignore list. Skipped')
            return
    if times == 0:
        DebugLog.info(f'Testing reaction "{rxnID}"')
    # finding demand for testing
    nextDemand = find_next_demand(model=model, check_rxn_id=rxnID, **kwargs)
    with suppress(ValueError):
        model.add_boundary(model.metabolites.get_by_id(nextDemand), "demand")
        model.reactions.get_by_id(f'DM_{nextDemand}').lower_bound = 0.2
        DebugLog.debug(f'Demand "DM_{nextDemand}" added')
    # Setting maximum times for recursion
    if times == len(model.reactions.get_by_id(rxnID).metabolites):
        msg = f'Reaction "{rxnID}" did not passed.'
        DebugLog.critical(msg)
        raise ValueError(msg)
    # answer must be reasonable and lie between given ranges
    # comparison must be using absolute values
    if not solution_range[0] <= abs(
            model.slim_optimize()) <= solution_range[1]:
        # Append to log
        DebugLog.debug(f'Reaction "{rxnID}" not in range')
        verify_sinks_for_rxn(
            model=model, rxnID=rxnID, **kwargs)
        # Recursive with 'extra' argument
        test_rxn_for_solution(
            model=model, rxnID=rxnID, solution_range=solution_range,
            times=times + 1, **kwargs)
    else:
        # if works, pass and return old objective
        DebugLog.debug(f'Reaction "{rxnID}" showed a feasible answer.')
        remove_boundary_if_not_model(
            model=model, metabolite=nextDemand, boundary="demand")
        verify_sinks_for_rxn(
            model=model, rxnID=rxnID, **kwargs)


def add_sequence(
        model: Model, sequence: list, avoid_list: Iterable = [], **kwargs):
    """From a sequence of Reaction objects, add reach Reaction to given model. It
    checks if new reaction does not break the metabolic system.

    :param model: Model to add the Reactions.
    :type model: Model
    :param sequence: List with Reaction objects
    :type sequence: list
    :raises TypeError: if Model is invalid
    :raises TypeError: if Reactions are not valid objects
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
            test_rxn_for_solution(model=model, rxnID=rxn.id, **kwargs)
            stopAndShowMassBalance(model=model, rxnID=rxn.id, **kwargs)
        else:
            # FIXME: avoid creating reaction
            DebugLog.warning(f'Reaction "{rxn.id}" was found in model')
    DebugLog.debug('Pathway added to Model')


def _add_graph_from_root(
        model: Model, root: Union[ET.Element, str], **kwargs):
    """For given root for a pathway, check, test and add all possible
    Reactions to given model.

    :param model: model to append pathway
    :type model: Model
    :param root: root or name of pathway
    :type root: Union[ET.Element, str]

    :Keyword Arguments:
        :param directory: Path to directory where data is located
        :type directory: Path
        :param database: Name for subdatabse, defaults to "META"
        :type database: str, optional
        :param comparment: location of the reactions
        :type comparment: compartment
        # TODO: Finish documentation and add more kwargs
        replacement_dict
        ignore_list
        avoid_list

    """
    # Retrieving and creating Pathway with Reactions
    # original_objective, original_direction = model.objective, \
    #     model.objective_direction
    graph = return_graph_from_root(root=root, **kwargs)
    for pathway in graph:
        sequence = list(create_reactions_for_iter(
            sequence=pathway, model=model, **kwargs))
        add_sequence(
            model=model, sequence=sequence, **kwargs)
    #TODO: Fix sink for metabolites in last sequence


def _add_graph_from_sequence(
        model: Model, sequence: Iterable, **kwargs):
        # TODO: add tests, and docstrings
    sequence = list(create_reactions_for_iter(
        model=model, sequence=sequence, **kwargs))
    add_sequence(
            model=model, sequence=sequence, **kwargs)


def add_graph_to_model(
        model: Model, graph: Union[list, str, ET.Element, set], **kwargs):
        # TODO: add tests, and docstrings
    if isinstance(graph, (ET.Element, str)):
        _add_graph_from_root(model=model, root=graph, **kwargs)
    elif isinstance(graph, (list, set)):
        _add_graph_from_sequence(model=model, sequence=graph, **kwargs)
