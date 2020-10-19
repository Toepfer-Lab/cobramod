#!/usr/bin/env python3
from typing import Union, Generator, Iterable, NamedTuple
from cobra import Model, Reaction, DictList
from cobramod.creation import check_mass_balance, _build_reaction
from cobramod.mod_parser import get_data
from cobramod.debug import debug_log
from itertools import chain
from contextlib import suppress
from pathlib import Path


class DataModel(NamedTuple):
    """
    A class to store old values of metabolic models
    """

    reactions: DictList
    metabolites: DictList
    boundary: DictList


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


def _replace_item(
    iterable: Iterable, replacement_dict: dict = {}
) -> Generator:
    """
    For an item in Iterable, replaces it for its corresponding value in
    given dictionary.

    Args:
        iterable (Iterable): sequence to modify
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.

    Yields:
        Generator: Either original keys or the replacements
    """
    for item in iterable:
        if item in set(chain.from_iterable(replacement_dict.keys())):
            debug_log.warning(
                f'Replacing "{item}" with "{replacement_dict[item]}"'
            )
            yield replacement_dict[item]
        else:
            yield item


def _remove_item(iterable: Iterable, avoid_list: Iterable = []) -> Generator:
    """
    Returns Generator of items that are not in the avoid list.
    """
    for item in iterable:
        if item in avoid_list:
            debug_log.warning(f'Avoiding root for "{item}"')
        else:
            yield item


def _return_verified_graph(
    vertex_dict: dict,
    vertex_set: set,
    avoid_list: Iterable = [],
    replacement_dict: dict = {},
) -> list:
    """
    Returns list with an ordered sequence for a dictionary of edges. If a
    vertex is missing in the sequence, they will be appended as a extra
    sequence.

    Args:
        vertex_dict (dict): dictionary with edges
        vertex_set (set): set with all members
        replacement_dict (dict, optional): original identifiers to be replaced.
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
    graph.insert(0, missing_edges)
    # Fixing sequences from graph
    graph = [
        list(
            _replace_item(iterable=sequence, replacement_dict=replacement_dict)
        )
        for sequence in graph
    ]
    graph = [
        list(_remove_item(iterable=sequence, avoid_list=avoid_list))
        for sequence in graph
    ]
    # Filtering NONES from .append()
    return list(filter(None, graph))


def return_graph_from_dict(
    data_dict: dict, avoid_list: Iterable = [], replacement_dict: dict = {}
) -> list:
    """
    Returns graph from given data. All sequences are in order.

    Args:
        data_dict (dict): pathway data saved in dictionary.
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.
        avoid_list (Iterable, Optional): original identifiers to avoid in
            graph. Defaults to [].
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.

    Returns:
        list: sequence with reaction identifiers in order
    """

    return _return_verified_graph(
        vertex_dict=data_dict["PATHWAY"],
        vertex_set=data_dict["SET"],
        avoid_list=avoid_list,
        replacement_dict=replacement_dict,
    )


def _create_reactions(
    sequence: Iterable,
    compartment: str,
    directory: Path,
    database: str,
    replacement_dict: dict = {},
    **kwargs,
) -> Generator:
    """
    For each identifier in the sequence, a Reaction will be created. It returns
    a generator.

    Deprecated form: _create_reactions_for_iter
    """

    debug_log.debug(f"Obtaining root for {sequence}")
    # From given list (which includes None values), retrieve only Reaction
    # Objects.
    for identifier in sequence:
        with suppress(KeyError):
            identifier = replacement_dict[identifier]
            debug_log.warning(
                f'Reaction "{identifier}" replaced by '
                f'"{replacement_dict[identifier]}".'
            )
        data_dict = get_data(
            identifier=identifier, database=database, directory=directory
        )
        yield _build_reaction(
            data_dict=data_dict,
            compartment=compartment,
            directory=directory,
            database=database,
            replacement_dict=replacement_dict,
            **kwargs,
        )


def _find_next_demand(
    model: Model, reaction_id: str, ignore_list: list = []
) -> str:
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
        Warning: if no metabolite was found

    Returns:
        str: metabolite identifier to create a demand reaction
    """
    tmp_rxn = model.reactions.get_by_id(reaction_id)
    # Checking reversibility
    # left --> right
    if model.reactions.get_by_id(reaction_id).upper_bound > abs(
        model.reactions.get_by_id(reaction_id).lower_bound
    ):
        tmp_list = [
            rxn.id for rxn in tmp_rxn.products if rxn.id not in ignore_list
        ]
    # left <-- right
    elif model.reactions.get_by_id(reaction_id).upper_bound < abs(
        model.reactions.get_by_id(reaction_id).lower_bound
    ):
        tmp_list = [
            rxn.id for rxn in tmp_rxn.reactants if rxn.id not in ignore_list
        ]
    # left <--> right
    elif model.reactions.get_by_id(reaction_id).upper_bound == abs(
        model.reactions.get_by_id(reaction_id).lower_bound
    ):
        # TODO: decide what to do (reversibility)
        # FIXME: isomers sometimes shows double demand
        tmp_list = [
            rxn.id for rxn in tmp_rxn.products if rxn.id not in ignore_list
        ]
    if len(tmp_list) == 0:
        # Nothing found
        raise Warning("No metabolite found to become a demand")
    else:
        debug_log.debug(
            f'Next demand selected for "{reaction_id}": "{tmp_list[0]}"'
        )
        return tmp_list[0]


def _has_demand(model: Model, metabolite: str) -> bool:
    """
    Returns True if model has a demand reaction for given metabolite
    identifier.
    """
    return any(
        [
            "DM_" in reaction.id
            for reaction in model.metabolites.get_by_id(metabolite).reactions
        ]
    )


def _remove_boundary_if_not_model(
    model: Model, metabolite: str, boundary: str
):
    """
    Removes given type of boundary reaction for a specified metabolite if found
    in the model.

    Args:
        model (Model): model to test
        metabolite (str): metabolite identifier in model
        boundary (str): type of boundary. Options: "sink", "demand"
    """
    type_prefix = {"sink": "SK_", "demand": "DM_"}
    if f"{type_prefix[boundary]}{metabolite}" in (
        reaction.id
        for reaction in model.metabolites.get_by_id(metabolite).reactions
    ):
        model.remove_reactions([f"{type_prefix[boundary]}{metabolite}"])
        debug_log.warning(
            f'{boundary.capitalize()} reaction for "{metabolite}" removed'
        )


def _less_equal_than_x_rxns(model: Model, metabolite: str, x: int) -> bool:
    """
    Returns True if given metabolite participates in less or equal than
    X reactions for given model.
    """
    reactions = model.metabolites.get_by_id(metabolite).reactions
    return len(reactions) <= x


def _fix_meta_for_boundaries(
    model: Model, metabolite: str, ignore_list: Iterable = []
):
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
    debug_log.debug(f'Checking "{metabolite}" in for sinks and demands')
    if metabolite in ignore_list:
        msg = f'Metabolite "{metabolite}" ignored'
        debug_log.warning(msg)
        raise Warning(msg)
    # Check for to add sinks
    if _has_demand(model=model, metabolite=metabolite):
        if _less_equal_than_x_rxns(model=model, metabolite=metabolite, x=2):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite), type="sink"
            )
            debug_log.warning(f'Sink reaction created for "{metabolite}"')
        else:
            _remove_boundary_if_not_model(
                model=model, metabolite=metabolite, boundary="sink"
            )
    else:
        if _less_equal_than_x_rxns(model=model, metabolite=metabolite, x=1):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite), type="sink"
            )
            debug_log.warning(f'Sink reaction created for "{metabolite}"')
        elif not _less_equal_than_x_rxns(
            model=model, metabolite=metabolite, x=2
        ):
            _remove_boundary_if_not_model(
                model=model, metabolite=metabolite, boundary="sink"
            )
    # remove demand if needed
    _remove_boundary_if_not_model(
        model=model, metabolite=metabolite, boundary="demand"
    )


def _verify_side_sinks_for_rxn(
    model: Model, rxn_id: str, side: str, ignore_list: Iterable = []
):
    """
    Checks for either the product or reactant side of a reactions, if
    participant-metabolites have enough sink reactions to produce a feasible
    answer. If necesary, it will creates or remove them.

    Args:
        model (Model): model to test.
        rxn_id (str): reaction identifier for given model.
        side (str): Side to verify. Options: "right", "left"
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []:

    Raises:
        TypeError: if ignore list is not iterable
        ValueError: if side is not in the options
    """
    if not isinstance(ignore_list, Iterable):
        raise TypeError("Ignore list is not iterable")
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
                model=model, metabolite=meta.id, ignore_list=ignore_list
            )


def verify_sinks_for_rxn(
    model: Model, rxn_id: str, ignore_list: Iterable = []
):
    """
    Verifies, creates or removes sink reactions for metabolites found in given
    reaction if certain conditions are met.

    .. hint:: If a metabolite has a demand reaction, it will create a sink\
    reaction if the total amount is equal or less than two. If no demand\
    reaction is found it will create a sink reaction if the total reactions is\
    1. Sink reactions will be delete if total reaction is either more than\
    equal 2, or 3.

    Args:
        model (Model): model to test
        rxn_id (str): reaction identifier for given model.
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []
    """
    # reactant side
    _verify_side_sinks_for_rxn(
        model=model, rxn_id=rxn_id, ignore_list=ignore_list, side="left"
    )
    # product side
    _verify_side_sinks_for_rxn(
        model=model, rxn_id=rxn_id, ignore_list=ignore_list, side="right"
    )


def _test_rxn_for_solution(
    model: Model,
    rxn_id: str,
    solution_range: tuple = (0.1, 1000),
    times: int = 0,
    ignore_list: list = [],
):
    """
    Checks if optimized objective function value for given model lies between
    a determinated range. Function is recursive and checks if sink reactions
    are enough or exceeded. It creates a demand reaction for reaction for the
    test and removes it, if necesary.

    Args:
        model (Model): model to test reaction
        rxn_id (str): reaction identifier in given model
        solution_range (tuple, optional): range of solution to pass the test.
            Defaults to (0.1, 1000).
        times (int, optional): Track of recursions. Defaults to 0.
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []

    Raises:
        ValueError: if solution is infeasible after many recursions. Depends
            from the amount of metabolites in the reaction
    """
    if rxn_id in ignore_list:
        debug_log.warning(f'Reaction "{rxn_id}" found in ignore list. Skipped')
        return
    # TODO: add counter
    if times == 0:
        debug_log.info(f'Testing reaction "{rxn_id}"')
    # finding demand for testing
    nextDemand = _find_next_demand(
        model=model, reaction_id=rxn_id, ignore_list=ignore_list
    )
    with suppress(ValueError):
        model.add_boundary(model.metabolites.get_by_id(nextDemand), "demand")
        model.reactions.get_by_id(f"DM_{nextDemand}").lower_bound = 0.2
        debug_log.debug(f'Demand "DM_{nextDemand}" added')
    # Setting maximum times for recursion
    if times == len(model.reactions.get_by_id(rxn_id).metabolites):
        msg = f'Reaction "{rxn_id}" did not passed.'
        debug_log.critical(msg)
        raise ValueError(msg)
    # answer must be reasonable and lie between given ranges
    # comparison must be using absolute values
    if (
        not solution_range[0]
        <= abs(model.slim_optimize())
        <= solution_range[1]
    ):
        # Append to log
        debug_log.debug(f'Reaction "{rxn_id}" not in range')
        verify_sinks_for_rxn(
            model=model, rxn_id=rxn_id, ignore_list=ignore_list
        )
        # Recursive with 'extra' argument
        _test_rxn_for_solution(
            model=model,
            rxn_id=rxn_id,
            solution_range=solution_range,
            times=times + 1,
            ignore_list=ignore_list,
        )
    else:
        # if works, pass and return old objective
        debug_log.debug(f'Reaction "{rxn_id}" showed a feasible answer.')
        _remove_boundary_if_not_model(
            model=model, metabolite=nextDemand, boundary="demand"
        )
        verify_sinks_for_rxn(
            model=model, rxn_id=rxn_id, ignore_list=ignore_list
        )


def _add_sequence(
    model: Model,
    sequence: list,
    avoid_list: Iterable = [],
    ignore_list: list = [],
    **kwargs,
):
    """
    From a sequence of Reaction objects, add each Reaction into given model. It
    checks if new reaction does not break the metabolic system.

    Args:
        model (Model): model to add reactions to.
        sequence (list): List with Reaction objects
        avoid_list (Iterable, optional): A sequence of formatted reactions to
            avoid adding to the model. This is usefull for long pathways,
            where X reactions need to be excluded. Defaults to [].
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to [].

    Keyword Arguments:
        solution_range (tuple, optional): range of solution to pass the test.
            Defaults to (0.1, 1000).

    Raises:
        TypeError: if reactions are not valid Reaction objects
    """
    # NOTE: temporal solution
    try:
        mass_kwargs = {
            "show_wrong": kwargs["show_wrong"],
            "stop_wrong": kwargs["stop_wrong"],
        }
    except KeyError:
        mass_kwargs = {}
    try:
        solution_range = kwargs["solution_range"]
    except KeyError:
        solution_range = (0.1, 1000)
    if not all((isinstance(reaction, Reaction) for reaction in sequence)):
        raise TypeError("Reactions are not valid objects. Check list")
    # only if not in model
    for rxn in sequence:
        if rxn.id in avoid_list:
            debug_log.warning(
                f'Reaction "{rxn.id}" found in avoid list. Skipping.'
            )
            continue
        if rxn.id not in model.reactions:
            model.add_reactions([rxn])
            debug_log.info(f'Reaction "{rxn.id}" added to model')
            _test_rxn_for_solution(
                model=model,
                rxn_id=rxn.id,
                ignore_list=ignore_list,
                solution_range=solution_range,
            )
            # This should be checked when creating the iter with Reaction
            check_mass_balance(model=model, rxn_id=rxn.id, **mass_kwargs)
        else:
            # FIXME: avoid creating reaction
            debug_log.info(f'Reaction "{rxn.id}" was found in model')
    debug_log.debug("Pathway added to Model")


def _from_data(
    model: Model,
    data_dict: dict,
    compartment: str,
    avoid_list: Iterable = [],
    replacement_dict: dict = {},
    ignore_list: list = [],
    **kwargs,
):
    """
    Adds root of a pathway into given model. Arguments explained above in
    `_add_sequence`.
    """
    # Retrieving and creating Pathway with Reactions
    # this uses replacment, avoid and xml
    graph = return_graph_from_dict(
        data_dict=data_dict,
        avoid_list=avoid_list,
        replacement_dict=replacement_dict,
    )
    for pathway in graph:
        sequence = list(
            _create_reactions(
                sequence=pathway,
                compartment=compartment,
                replacement_dict=replacement_dict,
                **kwargs,
            )
        )
        _add_sequence(
            model=model,
            sequence=sequence,
            ignore_list=ignore_list,
            avoid_list=avoid_list,
            **kwargs,
        )
        # TODO: Fix sink for metabolites in last sequence


def _from_sequence(
    model: Model,
    compartment: str,
    sequence: Iterable,
    database: str,
    directory: Path,
    avoid_list: Iterable = [],
    replacement_dict: dict = {},
    ignore_list: list = [],
    **kwargs,
):
    """
    Adds a sequence of identifiers to given model. Arguments explained
    above `_add_sequence`.
    """
    # TODO: add tests, and docstrings
    sequence = list(
        _create_reactions(
            sequence=sequence,
            compartment=compartment,
            database=database,
            replacement_dict=replacement_dict,
            directory=directory,
            **kwargs,
        )
    )
    _add_sequence(
        model=model,
        sequence=sequence,
        ignore_list=ignore_list,
        avoid_list=avoid_list,
        **kwargs,
    )


def _summary(old_values: DataModel, model: Model):
    print(
        f'{"-"*20}\n'
        f"Model: {model.id}\n"
        f"Original attributes:\n"
        f"Reactions: {len(old_values.reactions)}\n"
        f"Metabolites: {len(old_values.metabolites)}\n"
        f"Boundary reactions {len(old_values.boundary)}\n"
        f'{"-"*20}\n'
        f"New attributes:\n"
        f"Reactions: {len(model.reactions)}\n"
        f"Metabolites: {len(model.metabolites)}\n"
        f"Boundary reactions: {len(model.boundary)}\n"
        f'{"-"*20}\n'
    )
    # FIXME: ADD DIFF


def add_graph_to_model(
    model: Model,
    graph: Union[list, str, set],
    directory: Path,
    database: str,
    compartment: str,
    avoid_list: Iterable = [],
    replacement_dict: dict = {},
    ignore_list: list = [],
    **kwargs,
):
    """
    Adds given graph for a pathways or a sequence of reactions identifiers
    into given model. Data will be downloaded and structured according
    to the specified database.

    Args:
        model (Model): model to append graph
        graph (Union[list, str, set]): The identifier for the
            pathway or a iterator with the ids of the reactions.
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA", "KEGG"
        compartment (str): location of the reactions to take place. Defaults to
            cytosol "c"
        avoid_list (Iterable, optional): A sequence of formatted reactions to
            avoid adding to the model. This is usefull for long pathways,
            where X reactions need to be excluded.
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions.

    Keyword Arguments:
        show_wrong (bool): For each new reaction, if it is unbalance, it will
            show the difference for the metabolites. Defaults to True
        stop_wrong (bool): For each new reaction, if it is unbalance, it will
            stop and raise an error. Defaults to False
        solution_range (tuple, optional): range of solution to pass the test
            for each new reaction. Defaults to (0.1, 1000).

    Raises:
        TypeError: if model is invalid
    """
    if not isinstance(model, Model):
        raise TypeError("Model is invalid")
    old_values = DataModel(
        reactions=model.reactions.copy(),
        metabolites=model.metabolites.copy(),
        boundary=model.boundary.copy(),
    )
    if isinstance(graph, str):
        data_dict = get_data(
            directory=directory, identifier=graph, database=database
        )
        _from_data(
            model=model,
            data_dict=data_dict,
            compartment=compartment,
            avoid_list=avoid_list,
            replacement_dict=replacement_dict,
            ignore_list=ignore_list,
            directory=directory,
            database=database,
            **kwargs,
        )
        # FIXME: position
        _summary(old_values=old_values, model=model)
    elif isinstance(graph, (list, set)):
        _from_sequence(
            model=model,
            sequence=graph,
            compartment=compartment,
            avoid_list=avoid_list,
            directory=directory,
            database=database,
            replacement_dict=replacement_dict,
            ignore_list=ignore_list,
            **kwargs,
        )
        # FIXME: position
        _summary(old_values=old_values, model=model)
    else:
        raise TypeError("Argument graph must be iterable or a identifier")
