#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import logging
from typing import Union, Generator
from GLS import get_xml_from_biocyc, build_reaction_from_xml,\
    stopAndShowMassBalance
from cobra import Model, Reaction
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


def fix_single_reaction_dict(xmlRoot: ET.Element, rootDict: dict) -> dict:
    # These are single exceptions
    if len(rootDict) == 0:
        single_rxn = xmlRoot.find(
            "Pathway/reaction-list/Reaction").attrib["frameid"]
        rootDict = {single_rxn: single_rxn}
    return rootDict


def dictRxnForPathway(xmlRoot: Union[ET.Element, str], **kwargs) -> dict:
    """Returns a dictionary from the root of a .xml file, where the keys
    represent a Reaction and the value its predecessor in a pathway. i.e, the
    direction if from right to left

    :param xmlRoot: root of given pathway
    :type xmlRoot: Union[ET.Element, str]
    :raises TypeError: if root is not valid
    :return: dictionary with reactions and their predecessors
    :rtype: dict
    """
    rootDict = dict()
    if isinstance(xmlRoot, str):
        xmlRoot = get_xml_from_biocyc(bioID=xmlRoot, **kwargs)
    if not isinstance(xmlRoot, ET.Element):
        msg = f'xmlRoot is not a valid xml file.'
        DebugLog.error(msg)
        raise TypeError(msg)
    for rxnLine in xmlRoot.findall("*/reaction-ordering"):
        current = rxnLine.find("*/[@frameid]").attrib["frameid"]
        prior = rxnLine.find("predecessor-reactions/").attrib["frameid"]
        # If the direction of keys and values changes, then
        # many reactions would get lost. This way, even with
        # multiple compounds, pathways remain
        try:
            rootDict[current].append(prior)
        except KeyError:
            rootDict[current] = list(prior)
    name = xmlRoot.find("Pathway").attrib["frameid"]
    # NOTE: Add correspoding test
    rootDict = fix_single_reaction_dict(xmlRoot=xmlRoot, rootDict=rootDict)
    DebugLog.debug(f'Dictionary for pathway "{name}" succesfully created')
    return rootDict


def find_start_vertex(vertex_dict: dict) -> Generator:
    for key in vertex_dict.keys():
        if key not in [
            item for sublist in
                vertex_dict.values() for item in sublist]:
            yield key


def find_end_vertex(vertex_dict: dict) -> Generator:
    for value in [
            item for sublist in vertex_dict.values() for item in sublist]:
        if value not in [
            item for sublist in
                vertex_dict.keys() for item in sublist]:
            yield value


def fix_if_cycle(start: list, ends: list, vertex_dict: dict) -> tuple:
    if len(start) == 0:
        # TODO: check if randonmess affects model
        start = list(vertex_dict.keys())[0]
        if len(ends) == 0:
            ends = vertex_dict[start]
    return start, ends


def get_margin_vertex(vertex_dict: dict) -> tuple:
    """From a dictionary of Reactions IDs and their predecessors, returns the
    initial reaction and side/end-reactions.
    """
    vertex_start = list(
        find_start_vertex(vertex_dict=vertex_dict))
    vertex_ends = list(
        # could be multiple or cyclic
        find_end_vertex(vertex_dict=vertex_dict))
    # FIXME: cycles not working properly
    vertex_start, vertex_ends = fix_if_cycle(
        start=vertex_start, ends=vertex_ends, vertex_dict=vertex_dict)
    return vertex_start, vertex_ends


def give_path_graph(start_vertex: str, end_vertex_list: list, graph_dict: dict):
    # TO avoid cyclic
    # if start_vertex in end_vertex_list:
    #     end_vertex_list = []
    while start_vertex not in end_vertex_list:
        # try:
        next_vertex = graph_dict[start_vertex]
        for test in next_vertex:
            yield test
            # graph_dict.pop(test)
            start_vertex = test


def createPathway(
        xmlRoot: Union[ET.Element, str] = None, customDict: dict = None,
        **kwargs) -> list:
    """Creates from either a xmlRoot or custom dictionary a list with
    pathways. At least one argument should be passed.

    For custom dictionary direction is right to left


    :param xmlRoot: xmlroot of Pathway, defaults to None
    :type xmlRoot: Union[ET.Element, str], optional
    :param customDict: custom dictionary with specific entries, defaults
    to None
    :type customDict: dict, optional
    :raises TypeError: if no arguments are passed
    :raises Warning: if both arguments are passed. Only one arg is allowed
    :return: list with multiple list with Reactions IDs as strings
    :rtype: list
    """

    # creating dict
    if xmlRoot is None and customDict is None:
        raise TypeError(f'No arguments are passed')
    elif isinstance(xmlRoot, ET.Element) and isinstance(customDict, dict):
        raise Warning(f'Only one argument is accepted. Please re-run')
    elif customDict is None:
        rootDict = dictRxnForPathway(xmlRoot, **kwargs)  # completely normal
    elif isinstance(customDict, dict):  # custom dictionary
        rootDict = customDict
    # finding edges
    start, ends = get_edges_from_dict(rootDict)
    # Creating paths
    finalPaths = list()
    # TODO: make it more readable
    for path in ends:
        DebugLog.debug(f'Creating pathway ending with "{path}"')
        tmpList = list()  # sorted list
        tmpDict = rootDict.copy()  # copy is needed
        tmpList.append(path)  # end-metabolite
        # FIXME: Find a new solution to overcome cyclic graphs!!
        while True:
            try:
                if tmpDict[tmpList[-1]] in tmpDict.keys():
                    # When found, append to list and delete entry
                    # in dictionary.
                    tmpList.append(tmpDict[tmpList[-1]])
                    tmpDict.pop(tmpList[-2])
                # Must be the last metabolite
                elif tmpDict[tmpList[-1]] == start:
                    tmpList.append(tmpDict[tmpList[-1]])
                    break
            except KeyError:  # If only one key
                tmpList = ends
                break
        # Fixing index. Changing direction
        newSort = tmpList[::-1]
        DebugLog.debug(
            f'New Pathways is: '
            f'{newSort}')
        finalPaths.append(newSort)
    # else:
    #     # FIXME: This is an exceptio if only a reaction is in a pathway
    #     finalPaths.append([start])
    return finalPaths


def create_reactions_for_list(
        pathwayList: list, replaceNameDict: dict = None, **kwargs) -> list:
    """For each reaction ID in given list, it creates a Reaction object

    :param pathwayList: list with original reaction IDs
    :type pathwayList: list
    :param replaceNameDict: dictionary to replace names, defaults to None
    :type replaceNameDict: dict, optional
    :return: list with Reaction objects
    :rtype: list
    """
    search_kwargs = kwargs.copy()
    # TODO: write short FUN for this part
    for arg in ["model"]:
        try:
            del search_kwargs[arg]
        except KeyError:
            continue
    DebugLog.debug(f'Obtaining xmlRoot for {pathwayList}')
    # From given list (which includes None values), retrieve only Reaction
    # Objects.
    return [build_reaction_from_xml(
        xmlRoot=reaction, **kwargs) for reaction in pathwayList]


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
        model: Model, xmlRoot: Union[ET.Element, str], **kwargs):
    """For given xmlRoot for a pathway, check, test and add all possible
    Reactions to given model. 

    :param model: model to append pathway
    :type model: Model
    :param xmlRoot: root or name of pathway
    :type xmlRoot: Union[ET.Element, str]
    """	
    # Retrieving and creating Pathway with Reactions
    pathwayList = createPathway(xmlRoot=xmlRoot, **kwargs)
    for pathway in pathwayList:
        reactionList = create_reactions_for_list(
            pathwayList=pathway, model=model, **kwargs)
        addPathtoModel(model=model, listReactions=reactionList, **kwargs)
