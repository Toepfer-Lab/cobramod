#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import logging
from typing import Union, Generator, Iterable
from GLS import get_xml_from_biocyc, build_reaction_from_xml,\
    stopAndShowMassBalance
from cobra import Model, Reaction
from itertools import chain
# Creating corresponding Logs
# Format
DebugFormatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# Handler
DebugHandler = logging.FileHandler("debug.log", mode="a+")
DebugHandler.setFormatter(DebugFormatter)
# Log
DebugLog = logging.getLogger("DebugLog")
DebugLog.setLevel(logging.DEBUG)
# TODO!!: change name for proper python conventions


def fix_single_reaction_dict(root: ET.Element, rootDict: dict) -> dict:
    # TODO: verify
    # These are single exceptions
    if len(rootDict) == 0:
        single_rxn = root.find(
            "Pathway/reaction-list/Reaction").attrib["frameid"]
        rootDict = {single_rxn: single_rxn}
    return rootDict


def create_vertex_dict(root: Union[ET.Element, str], **kwargs) -> dict:
    """Returns a dictionary from the root of a XML file, where the keys
    represent a Reaction and the value its predecessor in a pathway. i.e, the
    direction if from right to left

    :param root: root of given pathway
    :type root: Union[ET.Element, str]
    :raises TypeError: if root is not valid
    :return: dictionary with reactions and their predecessors
    :rtype: dict
    """
    rootDict = dict()
    if isinstance(root, str):
        root = get_xml_from_biocyc(bioID=root, **kwargs)
    if not isinstance(root, ET.Element):
        msg = f'root is not a valid xml file.'
        DebugLog.error(msg)
        raise TypeError(msg)
    for rxnLine in root.findall("*/reaction-ordering"):
        current = rxnLine.find("*/[@frameid]").attrib["frameid"]
        prior = rxnLine.find("predecessor-reactions/").attrib["frameid"]
        # If the direction of keys and values changes, then
        # many reactions would get lost. This way, even with
        # multiple compounds, pathways remain
        # NOTE: check for behaviour
        # Replacing values produces cuts with are needed to avoid cyclic
        # No reactions are missed
        rootDict[current] = prior
    name = root.find("Pathway").attrib["frameid"]
    # NOTE: Add correspoding test
    rootDict = fix_single_reaction_dict(root=root, rootDict=rootDict)
    DebugLog.debug(f'Dictionary for pathway "{name}" succesfully created')
    return rootDict


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
    end_vertex = verify_return_end_vertex(vertex_dict=vertex_dict)
    start_vertex = find_start_vertex(vertex_dict=vertex_dict)
    # a generator with list, otherwise items get lost
    return (list(build_graph(
        start_vertex=start, vertex_dict=vertex_dict)) for start in start_vertex)


def return_missing_edges(complete_edges:Iterable, graph:Iterable) -> Generator:
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


def return_verified_graph(vertex_dict:dict) -> list:
    """Returns list with ordered sequences for a edge dictionary. If missing
    vertex are missing in the sequences, they will be appended in the list

    :param vertex_dict: dictionary of edges
    :type vertex_dict: dict
    :return: verified list with sequences
    :rtype: list
    """
    # Backup for comparison
    vertex_dict_copy = vertex_dict.copy()
    all_edges = set(
            chain.from_iterable(vertex_dict_copy.items()))
    # TODO: check for single reaction
    graph = list(get_graph(vertex_dict=vertex_dict))
    missing_edges = list(return_missing_edges(
        graph=set(chain.from_iterable(graph)),
        complete_edges=all_edges))
    graph.append(missing_edges)
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
        raise AttributeError(f'Given XML root is invalid.')
    vertex_dict = create_vertex_dict(root=root, **kwargs)
    return return_verified_graph(vertex_dict=vertex_dict)


def create_reactions_for_iter(
        pathwayList: list, replaceNameDict: dict = None,
        **kwargs) -> Generator:
    """For each reaction ID in given list, it creates a Reaction object

    :param pathwayList: list with original reaction IDs
    :type pathwayList: list
    :param replaceNameDict: dictionary to replace names, defaults to None
    :type replaceNameDict: dict, optional
    :return: list with Reaction objects
    :rtype: list
    """
    search_kwargs = kwargs.copy()
    # TODO: verify is this needed
    for arg in ["model"]:
        try:
            del search_kwargs[arg]
        except KeyError:
            continue
    DebugLog.debug(f'Obtaining root for {pathwayList}')
    # From given list (which includes None values), retrieve only Reaction
    # Objects.
    return (build_reaction_from_xml(
        root=reaction, **kwargs) for reaction in pathwayList)


def find_next_demand(
        model: Model, toCheckRxnID: str, ignoreList: list = [],
        **kwargs) -> str:
    """Returns first metabolites found either in the product or reactant side
    of given reaction. Reversibility of the reaction is taken into
    consideration. A list with metabolites IDs can be passed to be ignored

    :param model: Model where the reaction is located
    :type model: Model
    :param toCheckRxnID: ID of the reaction in the model
    :type toCheckRxnID: str
    :param ignoreList: list with metabolites IDs to ignore, defaults to []
    :type ignoreList: list, optional
    :raises TypeError: if Model is not valid
    :raises Warning: if no metabolite was found
    :return: ID of next metabolite to become demand
    :rtype: str
    """
    if not isinstance(model, Model):
        raise TypeError(f'Model is invalid')
    tmpRxn = model.reactions.get_by_id(toCheckRxnID)
    # left --> right
    if model.reactions.get_by_id(toCheckRxnID).upper_bound > abs(
            model.reactions.get_by_id(toCheckRxnID).lower_bound):
        tmpList = [
            rxn.id for rxn in tmpRxn.products if rxn.id not in ignoreList]
    # left <-- right
    elif model.reactions.get_by_id(toCheckRxnID).upper_bound < abs(
            model.reactions.get_by_id(toCheckRxnID).lower_bound):
        tmpList = [
            rxn.id for rxn in tmpRxn.reactants if rxn.id not in ignoreList]
    # left <--> right
    elif model.reactions.get_by_id(toCheckRxnID).upper_bound == abs(
            model.reactions.get_by_id(toCheckRxnID).lower_bound):
        # TODO: decide what to do (reversibility)
        tmpList = [
            rxn.id for rxn in tmpRxn.products if rxn.id not in ignoreList]
    if len(tmpList) == 0:
        # Nothing found
        raise Warning(
            f'No metabolite found to become a demand')
    else:
        DebugLog.info(
            f'Next demand selected for "{toCheckRxnID}": "{tmpList[0]}"')
        return tmpList[0]


def isSink(model: Model, rxnID: str) -> bool:
    """Returns True if given reaction ID has a sink reaction.

    :param model: model to check for sinks
    :type model: Model
    :param rxnID: ID of reaction in given Model
    :type rxnID: str
    :return: True if found, False is not
    :rtype: bool
    """
    allSinks = [reaction.id for reaction in model.sinks]
    return rxnID in allSinks


def has_demand(model: Model, rxnID: str) -> bool:
    """Returns True if given reaction ID has a sink reaction.

    :param model: model to check for sinks
    :type model: Model
    :param rxnID: ID of reaction in given Model
    :type rxnID: str
    :return: True if found, False is not
    :rtype: bool
    """
    allSinks = [reaction.id for reaction in model.sinks]
    return rxnID in allSinks


def need_demand_metabolite(
        model: Model, metaID: str) -> bool:
    return not len([
        reaction.id for reaction in model.metabolites.get_by_id(
            metaID).reactions]) > 2


def remove_boundary_if_not_model(
        model: Model, metabolite: str, boundary: str = "sink"):
    typeDict = {"sink": "SK_", "demand": "DM_"}
    if f'{typeDict[boundary]}{metabolite}' in [
        reaction.id for reaction in model.metabolites.get_by_id(
            metabolite).reactions]:
        model.remove_reactions([f'{typeDict[boundary]}{metabolite}'])
        DebugLog.warning(
            f'{boundary.capitalize()} reaction for "{metabolite}" removed')


def more_than_two_reaction(model: Model, metaID: str) -> bool:
    """Returns if reactions where given metabolite participates, is larger than
    bool""" #FIXME:
    reactions = model.metabolites.get_by_id(metaID).reactions
    return len(reactions) > 2


def less_than_two_reaction(model: Model, metaID: str) -> bool:
    #FIXME:
    reactions = model.metabolites.get_by_id(metaID).reactions
    return len(reactions) < 2


def check_if_boundary(model: Model, metabolite: str) -> bool:
    list_no_exchanges = [demand.id for demand in model.demands]\
        + [sink.id for sink in model.sinks]
    return any([
        f'DM_{metabolite}' in list_no_exchanges,
        f'SK_{metabolite}' in list_no_exchanges])


def fix_meta_for_boundaries(
        model: Model, metabolite: str, ignoreList: list = [], **kwargs):
    DebugLog.debug(f'Checking "{metabolite}" in for sinks and demands')
    if metabolite in ignoreList:
        msg = f'Metabolite "{metabolite}" ignored'
        DebugLog.warning(msg)
        raise Warning(msg)
    # check if boundaries besides exchange
    if check_if_boundary(model=model, metabolite=metabolite):
        if more_than_two_reaction(model=model, metaID=metabolite):
            remove_boundary_if_not_model(
                model=model, metabolite=metabolite)
        if not need_demand_metabolite(model=model, metaID=metabolite):
            remove_boundary_if_not_model(
                model=model, metabolite=metabolite, boundary="demand")
    else:
        if less_than_two_reaction(model=model, metaID=metabolite):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite),
                type="sink",)
            DebugLog.warning(f'Sink reaction created for "{metabolite}"')


def tryCreateOrRemoveSinkForSides(
        model: Model, rxnID: str, side: str = "right", ignoreList: list = []):
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
    :param ignoreList: list with metabolites to be ignored, defaults to []
    :type ignoreList: list, optional
    :raises TypeError: if model is invalid
    :raises TypeError: if ignoreList is not passed as a list
    :raises ValueError: If 'side' is not properly passed
    """
    if not isinstance(ignoreList, list):
        raise TypeError(f'IgnoreList is not passed as a list')
    if not isinstance(model, Model):
        raise TypeError(f'Model is invalid')
    if side == "right":
        metabolites = model.reactions.get_by_id(rxnID).products
    elif side == "left":
        metabolites = model.reactions.get_by_id(rxnID).reactants
    else:
        raise ValueError(f'Only valid options are "right" and "left"')
    for meta in metabolites:
        # create or remove
        try:
            fix_meta_for_boundaries(
                model=model, metabolite=meta.id, ignoreList=ignoreList)
        # For ignored metabolites
        except Warning:
            pass


# TODO: check if kwargs are needed
def createAndCheckSinksForRxn(
        model: Model, rxnID: str, ignoreList: list = [], **kwargs):
    """Checks, creates or remove sink reactions for metabolites in given
    Reaction

    :param model: model with Reactions
    :type model: Model
    :param rxnID: ID of the reaction in model
    :type rxnID: str
    :param ignoreList: metabolites to be ignored, defaults to []
    :type ignoreList: list, optional
    """
    # reactant side
    tryCreateOrRemoveSinkForSides(
        model=model, rxnID=rxnID, ignoreList=ignoreList, side="left")
    # product side
    tryCreateOrRemoveSinkForSides(
        model=model, rxnID=rxnID, ignoreList=ignoreList, side="right")
 

def testReactionForSolution(
        model: Model, rxnID: str, solutionRange: tuple = (0.1, 1000),
        times: int = 0, **kwargs):
    """Checks if optimized objective function value of given model lies between
    a determinated range. Function is recursive and checks if sink reactions
    are enough or exceded. It creates a demand reaction for reaction for the
    test and removes it, if necesary.

    :param model: model to test
    :type model: Model
    :param rxnID: reaction ID to test
    :type rxnID: str
    :param solutionRange: range of solution, defaults to (0.1, 1000)
    :type solutionRange: tuple, optional
    :param times: For recursion, displays how many times the FUN called itself, defaults to 0
    :type times: int, optional
    :raises ValueError: if solution is infeasible after many recursions
    """
    # TODO: add errors for wrong MODEL or rxnid
    if times == 0:
        DebugLog.debug(f'Testing reaction "{rxnID}"')
    model.objective = rxnID
    # finding demand for testing
    nextDemand = find_next_demand(model=model, toCheckRxnID=rxnID, **kwargs)
    # New objective.
    try:
        model.add_boundary(model.metabolites.get_by_id(nextDemand), "demand")
        model.reactions.get_by_id(f'DM_{nextDemand}').lower_bound = 0.1
        DebugLog.debug(f'Demand "DM_{nextDemand}" added')
    except ValueError:
        # boundary already in Model
        pass
    # changing direction for negative solutions
    model.objective_direction = {
        "direction": "min" if abs(
            model.reactions.get_by_id(rxnID).lower_bound) > abs(
                model.reactions.get_by_id(rxnID).upper_bound) else "max"}.get(
                    "direction")
    if times == len(model.reactions.get_by_id(rxnID).metabolites):
        msg = f'Reaction "{rxnID}" did not passed.'
        DebugLog.critical(msg)
        raise ValueError(msg)  # check for errors
    # answer must be reasonable and lie between given ranges
    # comparison must be using absolute values
    if not solutionRange[0] <= abs(model.slim_optimize()) <= solutionRange[1]:
        # Append to log
        DebugLog.debug(f'Reaction "{rxnID}" not in range')
        createAndCheckSinksForRxn(
            model=model, rxnID=rxnID, **kwargs)
        # Recursive with 'extra' argument
        testReactionForSolution(
            model=model, rxnID=rxnID, solutionRange=solutionRange,
            times=times + 1, **kwargs)
    else:
        # if works, pass and return old objective
        DebugLog.debug(f'Reaction "{rxnID}" showed a feasible answer.')
        # for demand
        # FIXME: a metabolite might get new reactions. Recheck for demand
        if more_than_two_reaction(model=model, metaID=nextDemand):
            model.remove_reactions([f'DM_{nextDemand}'])
            DebugLog.warning(f'Demand "DM_{nextDemand}" removed')
        createAndCheckSinksForRxn(
            model=model, rxnID=rxnID, **kwargs)
        # else:
        #     pass


def addPathtoModel(model: Model, listReactions: list, **kwargs):
    """From a list of Reaction objects, add reach Reaction to given model. It
    checks if new reaction does not break the metabolic system.

    :param model: Model to add the Reactions.
    :type model: Model
    :param listReactions: List with Reaction objects
    :type listReactions: list
    :raises TypeError: if Model is invalid
    :raises TypeError: if Reactions are not valid objects
    """
    if not isinstance(model, Model):
        raise TypeError(f'Model is invalid')
    if not all([isinstance(reaction, Reaction) for reaction in listReactions]):
        raise TypeError(f'Reactions are not valid objects. Check list')
    # only if not in model
    for rxn in listReactions:
        if rxn not in model.reactions:
            model.add_reactions([rxn])
            DebugLog.info(f'Reaction "{rxn.id}" added to model')
            testReactionForSolution(model=model, rxnID=rxn.id, **kwargs)
            stopAndShowMassBalance(model=model, rxnID=rxn.id, **kwargs)
    DebugLog.debug(f'Pathway added to Model')


def testAndAddCompletePathway(
        model: Model, root: Union[ET.Element, str], **kwargs):
    """For given root for a pathway, check, test and add all possible
    Reactions to given model. 

    :param model: model to append pathway
    :type model: Model
    :param root: root or name of pathway
    :type root: Union[ET.Element, str]
    """	
    # Retrieving and creating Pathway with Reactions
    pathwayList = createPathway(root=root, **kwargs)
    for pathway in pathwayList:
        reactionList = create_reactions_for_iter(
            pathwayList=pathway, model=model, **kwargs)
        addPathtoModel(model=model, listReactions=reactionList, **kwargs)
