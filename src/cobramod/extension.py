#!/usr/bin/env python3
"""Pathway extension

This module handles the addition of reactions as a pathway in to
a model and the corresponding test that come with it.

Most important functions:
- add_pathway: Adds a pathway or multiple reactions into a model.
- test_result: Checks that given reaction in a model is active and gives a
flux.
"""
from contextlib import suppress
from pathlib import Path
from typing import Union, Generator, Iterable

from cobra import Model, Reaction
from requests import HTTPError

from cobramod.creation import create_object
from cobramod.debug import debug_log
from cobramod.error import NotInRangeError
from cobramod.graph import return_graph_from_dict
from cobramod.mod_parser import get_data
from cobramod.utils import get_DataList, get_basic_info, check_to_write
from cobramod.pathway import Pathway


def _create_reactions(
    sequence: list,
    compartment: str,
    directory: Path,
    database: str,
    show_imbalance: bool,
    stop_imbalance: bool,
    replacement: dict,
    model: Model,
    model_id: str,
) -> Generator:
    """
    For each identifier in the sequence, a Reaction will be created. It returns
    a generator.

    .. hint:: Hyphens will become underscores. Double hyphens become single\
    underscores.

    Args:
        sequence (list): Identifiers for the reactions.
        directory (Path): Path to directory where data is located.
        database (str): Name of the database. Options are: "META", "ARA",
            "KEGG", "BIGG"
        compartment (str): Location of the reaction.
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        model_id (str): Exclusive for BIGG. Retrieve object from specified
            model.

    Returns:
        Generator: A generator to that gives the reactions.
    """

    debug_log.debug(f"Obtaining data for {sequence}")
    # From given list (which could include None values), retrieve only Reaction
    # Objects.
    for identifier in sequence:
        with suppress(KeyError):
            identifier = replacement[identifier]
            debug_log.warning(
                f'Reaction "{identifier}" replaced by '
                f'"{replacement[identifier]}".'
            )
        # Check if translation is available
        yield create_object(
            identifier=identifier,
            directory=directory,
            database=database,
            compartment=compartment,
            replacement=replacement,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            model=model,
            model_id=model_id,
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
    tmp_reaction = model.reactions.get_by_id(reaction_id)
    # FIXME: find a better approach
    # Checking reversibility
    # left --> right
    if model.reactions.get_by_id(reaction_id).upper_bound > abs(
        model.reactions.get_by_id(reaction_id).lower_bound
    ):
        tmp_list = [
            reaction.id
            for reaction in tmp_reaction.products
            if reaction.id not in ignore_list
        ]
    # left <-- right
    elif model.reactions.get_by_id(reaction_id).upper_bound < abs(
        model.reactions.get_by_id(reaction_id).lower_bound
    ):
        tmp_list = [
            reaction.id
            for reaction in tmp_reaction.reactants
            if reaction.id not in ignore_list
        ]
    # left <--> right
    elif model.reactions.get_by_id(reaction_id).upper_bound == abs(
        model.reactions.get_by_id(reaction_id).lower_bound
    ):
        # TODO: decide what to do (reversibility)
        # FIXME: isomers sometimes shows double demand
        tmp_list = [
            reaction.id
            for reaction in tmp_reaction.products
            if reaction.id not in ignore_list
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


def _remove_boundary(model: Model, metabolite: str, boundary: str):
    """
    Removes given type of boundary reaction for a specified metabolite if it
    is found in the model.

    Args:
        model (Model): model to examine
        metabolite (str): metabolite identifier in model to examine.
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


def _less_equal_than_x_reactions(
    model: Model, metabolite: str, x: int
) -> bool:
    """
    Returns True if given metabolite participates in less or equal than
    X reactions for given model.
    """
    reactions = model.metabolites.get_by_id(metabolite).reactions
    return len(reactions) <= x


def _verify_boundary(model: Model, metabolite: str, ignore_list: Iterable):
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
    debug_log.debug(f'Checking "{metabolite}" for sinks and demands')
    if metabolite in ignore_list:
        msg = f'Metabolite "{metabolite}" ignored'
        debug_log.warning(msg)
        raise Warning(msg)
    # Check for to add sinks. If demand reaction is found and metabolite has
    # less than two reactions, then create sink, otherwise remove any sink.
    if _has_demand(model=model, metabolite=metabolite):
        if _less_equal_than_x_reactions(
            model=model, metabolite=metabolite, x=2
        ):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite), type="sink"
            )
            debug_log.warning(f'Sink reaction created for "{metabolite}"')
        else:
            _remove_boundary(
                model=model, metabolite=metabolite, boundary="sink"
            )
    # Otherwise, if it has no demand, add the sink if it has only one reaction
    # and remove sink if more than 2 reactions.
    # TODO: check for else behaviour
    else:
        if _less_equal_than_x_reactions(
            model=model, metabolite=metabolite, x=1
        ):
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite), type="sink"
            )
            debug_log.warning(f'Sink reaction created for "{metabolite}"')
        elif not _less_equal_than_x_reactions(
            model=model, metabolite=metabolite, x=2
        ):
            _remove_boundary(
                model=model, metabolite=metabolite, boundary="sink"
            )
    # Remove demand if needed. This part is always necessary since the extra
    # demands would add extra flux into the system.
    _remove_boundary(model=model, metabolite=metabolite, boundary="demand")


def _fix_side(
    model: Model, reaction: str, side: str, ignore_list: Iterable = []
):
    """
    Checks for either the product or reactant side of a reactions, if
    participant-metabolites have enough sink reactions to produce a feasible
    answer. If necesary, it will creates or remove them.

    Args:
        model (Model): model to test.
        reaction_id (str): reaction identifier for given model.
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
        metabolites = model.reactions.get_by_id(reaction).products
    elif side == "left":
        metabolites = model.reactions.get_by_id(reaction).reactants
    else:
        raise ValueError('Only valid options are "right" and "left"')
    for meta in metabolites:
        # create or remove
        with suppress(Warning):
            _verify_boundary(
                model=model, metabolite=meta.id, ignore_list=ignore_list
            )


def _verify_sinks(model: Model, reaction: str, ignore_list: Iterable):
    """
    Verifies, creates or removes sink reactions for metabolites found in given
    reaction if certain conditions are met.

    .. hint:: If a metabolite has a demand reaction, it will create a sink\
    reaction if the total amount is equal or less than two. If no demand\
    reaction is found it will create a sink reaction if the total reactions is\
    1. Sink reactions will be delete if total reaction is either more than\
    equal 2, or 3.

    Args:
        model (Model): model to examine.
        reaction (str): reaction identifier for given model.
        ignore_list (Iterable): A sequence of formatted metabolites
            to ignore when testing new reactions.
    """
    # reactant side
    _fix_side(
        model=model, reaction=reaction, ignore_list=ignore_list, side="left"
    )
    # product side
    _fix_side(
        model=model, reaction=reaction, ignore_list=ignore_list, side="right"
    )


def test_result(
    model: Model,
    reaction: str,
    minimum: float = 0.1,
    times: int = 0,
    ignore_list: list = [],
):
    """
    Checks if optimized objective function value for given model has a minimum
    value. Function is recursive and checks if sink reactions
    are enough or exceeded. It creates a demand reaction for reaction for the
    test and removes it, if necesary.

    Args:
        model (Model): model, where reaction is located
        reaction (str): reaction identifier in model to test.
        minium (float, optional): Minium value for test to pass.
            Defaults to 0.1 .
        times (int, optional): Track of recursions. Defaults to 0.
        ignore_list (Iterable, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []

    Raises:
        NotInRangeError: if solution is infeasible after many recursions.
            Depends from the amount of metabolites in the reaction.
    """
    if reaction in ignore_list:
        debug_log.warning(
            f'Reaction "{reaction}" found in ignore list. Skipped'
        )
        return
    # TODO: add counter
    if times == 0:
        debug_log.info(f'Testing reaction "{reaction}"')
    # finding demand for testing
    next_demand = _find_next_demand(
        model=model, reaction_id=reaction, ignore_list=ignore_list
    )
    # The first step is to add a demand for check that the reaction has an
    # active flux.
    with suppress(ValueError):
        model.add_boundary(model.metabolites.get_by_id(next_demand), "demand")
        model.reactions.get_by_id(f"DM_{next_demand}").lower_bound = 0.1
        debug_log.debug(f'Demand "DM_{next_demand}" added')
    # Setting maximum times for recursion
    if times == len(model.reactions.get_by_id(reaction).metabolites):
        raise NotInRangeError(reaction=reaction)
    # answer must be reasonable
    # comparison must be using absolute values
    if not minimum <= abs(model.slim_optimize()):
        # Append to log
        debug_log.debug(f'Reaction "{reaction}" not in range')
        _verify_sinks(model=model, reaction=reaction, ignore_list=ignore_list)
        # Recursive with 'extra' argument
        test_result(
            model=model,
            reaction=reaction,
            minimum=minimum,
            times=times + 1,
            ignore_list=ignore_list,
        )
    else:
        # if works, pass and return old objective
        debug_log.debug(f'Reaction "{reaction}" showed a feasible answer.')
        # Remove old demand
        _remove_boundary(
            model=model, metabolite=next_demand, boundary="demand"
        )
        _verify_sinks(model=model, reaction=reaction, ignore_list=ignore_list)


def _add_sequence(
    model: Model,
    identifier: str,
    sequence: list,
    avoid_list: list,
    ignore_list: list,
    minimum: float,
):
    """
    From a sequence of Reaction objects, add each Reaction into given model. It
    checks if new reactions do not break the optimized value. All reactions
    will be added to a common COBRApy Group.

    Args:
        model (Model): Model to expand.
        identifier (str): Common :func:`cobramod.pathway.Pathway` identifier.
        sequence (list): List with :func:`cobra.core.reaction.Reaction` objects
        avoid_list (list): A sequence of formatted reactions to
            avoid adding to the model. This is useful for long pathways,
            where X reactions need to be excluded.
        ignore_list (list): A sequence of formatted metabolites to skip when
            testing, and/or reactions that should be added but not tested.
            This is useful for long cyclical pathways.
        minimun (float): Minimum optimized value to pass in every single test.

    Raises:
        TypeError: if reactions are not valid Reaction objects
    """
    if not all((isinstance(reaction, Reaction) for reaction in sequence)):
        raise TypeError("Reactions are not valid objects. Check list")
    # Either create a Pathway or obtain the correct Pathway.
    pathway = Pathway(id=identifier)
    if identifier in [group.id for group in model.groups]:
        pathway = model.groups.get_by_id(identifier)
    # Add sequence to model
    for reaction in sequence:
        # This reaction will not be added.
        if reaction.id in avoid_list:
            debug_log.warning(
                f'Reaction "{reaction.id}" found in avoid list. Skipping.'
            )
            continue
        # FIXME: avoid if statements for perfomance improvement
        if reaction.id not in [reaction.id for reaction in model.reactions]:
            model.add_reactions([reaction])
            debug_log.info(f'Reaction "{reaction.id}" added to model')
            debug_log.debug(
                f'Reaction "{reaction.id}" added to group "{pathway.id}"'
            )
            if reaction.id not in ignore_list:
                test_result(
                    model=model,
                    reaction=reaction.id,
                    ignore_list=ignore_list,
                    minimum=minimum,
                )
            else:
                debug_log.warning(f'Test for reaction "{reaction.id}" skipped')
            # Add to pathway only if reaction was not previously in the model.
            pathway.add_members(new_members=[reaction])
        else:
            # FIXME: avoid creating reaction
            debug_log.info(
                f'Reaction "{reaction.id}" was found in model. Skipping.'
            )
    # TODO: Add space
    debug_log.debug(f'Reactions added to group "{pathway.id}"')
    if identifier not in [group.id for group in model.groups]:
        model.add_groups(group_list=[pathway])
        debug_log.debug("Pathway added to Model")


def _from_data(
    model: Model,
    data_dict: dict,
    directory: Path,
    database: str,
    compartment: str,
    avoid_list: list,
    replacement: dict,
    ignore_list: list,
    minimum: float,
    stop_imbalance: bool,
    show_imbalance: bool,
    model_id: str,
):
    """"
    Adds a pathway into given model from a dictinary with the information of
    the pathway into given model from a dictinary with the information of
    the model.

    Args:
        model (Model): Model to expand.
        data_dict (dict): Dictinary with the information for the pathway.
        directory (Path): Path for directory to stored and retrieve data.
        database (str): Name of the database to search for reactions and
            metabolites.
        compartment: Location of the reactions.

    Arguments for complex pathways:
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
        avoid_list (Iterable): A sequence of formatted reactions to
            avoid adding to the model. This is useful for long pathways,
            where X reactions need to be excluded.
        ignore_list (list): A sequence of formatted metabolites to skip when
            testing, and/or reactions that should be added but not tested.
            This is useful for long cyclical pathways.
        minimun (float): Minimum optimized value to pass in every single test.

    Arguments for utilities:
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathway are not available.
            Defaults to: "universal"
    """
    # A graph can have multiple routes, depending on how many end-metabolites.
    graph = return_graph_from_dict(
        data_dict=data_dict, avoid_list=avoid_list, replacement=replacement
    )
    for sequence in graph:
        # Data storage is handled by method
        sequence = list(
            _create_reactions(
                sequence=sequence,
                compartment=compartment,
                directory=directory,
                database=database,
                replacement=replacement,
                stop_imbalance=stop_imbalance,
                show_imbalance=show_imbalance,
                model=model,
                model_id=model_id,
            )
        )
        # Add to model
        _add_sequence(
            model=model,
            identifier=data_dict["ENTRY"],
            sequence=sequence,
            ignore_list=ignore_list,
            avoid_list=avoid_list,
            minimum=minimum,
        )
        # TODO: Fix sink for metabolites in last sequence


def _from_sequence(
    model: Model,
    identifier: str,
    sequence: list,
    database: str,
    compartment: str,
    directory: Path,
    avoid_list: list,
    replacement: dict,
    ignore_list: list,
    stop_imbalance: bool,
    show_imbalance: bool,
    minimum: float,
    model_id: str,
):
    """
    Adds a sequence of identifiers to given model. It will automatically test
    for valid optimized values.

    Args:
        model (Model): Model to expand.
        identifier (str): Common :func:`cobramod.pathway.Pathway` identifier.
        sequence (list): List reaction identifiers.
        database (str): Name of the database.
        compartment: Location of the reactions.
        directory (Path): Path for directory to stored and retrieve data.

    Arguments for complex pathways:
        avoid_list (list): A sequence of formatted reactions to
            avoid adding to the model. This is useful for long pathways,
            where X reactions need to be excluded.
        ignore_list (list): A sequence of formatted metabolites to skip when
            testing, and/or reactions that should be added but not tested.
            This is useful for long cyclical pathways.
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
        minimun (float): Minimum optimized value to pass in every single test.

    Arguments for utilities:
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathway are not available.
    """
    # Data storage is managed by the function
    sequence = list(
        _create_reactions(
            sequence=sequence,
            compartment=compartment,
            database=database,
            replacement=replacement,
            directory=directory,
            stop_imbalance=stop_imbalance,
            show_imbalance=show_imbalance,
            model=model,
            model_id=model_id,
        )
    )
    # Append to model
    _add_sequence(
        model=model,
        sequence=sequence,
        identifier=identifier,
        ignore_list=ignore_list,
        avoid_list=avoid_list,
        minimum=minimum,
    )


def add_pathway(
    model: Model,
    pathway: Union[list, str],
    directory: Path,
    database: str,
    compartment: str,
    group: str = "custom_group",
    avoid_list: list = [],
    replacement: dict = {},
    ignore_list: list = [],
    filename: Path = None,
    summary: bool = True,
    minimum: float = 0.1,
    stop_imbalance: bool = False,
    show_imbalance: bool = True,
    model_id: str = "universal",
):
    """
    Adds given graph for a pathway identifier or a sequence of reactions
    identifiers into given model. Data will be downloaded and structured
    according to the specified database.

    Args:
        model (Model): Model to expand.
        pathway (list, str): Sequence of reaction identifiers or
            identifier for a pathway. Examples: ["RXN-2206", "RXN-207"] or
            "PWY-886"
        directory (Path): Path for directory to stored and retrieve data.
        database (str): Name of the database.
        compartment: Location of the reactions.
        group (str, optional): Common :func:`cobramod.pathway.Pathway`
            identifier. Defaults to "custom_group"

    Arguments for complex pathways:
        avoid_list (list, optional): A sequence of formatted reactions to avoid
            adding to the model. This is useful for long pathways, where X
            reactions need to be excluded.
        replacement (dict, optional): Original identifiers to be replaced.
            Values are the new identifiers.
        ignore_list (list): A sequence of formatted metabolites to skip when
            testing, and/or reactions that should be added but not tested.
            This is useful for long cyclical pathways.
        minimum (float, optional): Minimum optimized value to pass in every
            single test. Defaults to 0.1

    Arguments for utilities:
        filename (Path): Location for the summary. Defaults to "summary.txt" in
            the current working directory.
        summary (bool, optional): True to write summary in file. Defaults to
            True.
        stop_imbalance (bool, optional): If unbalanced reaction is found, stop
            process. Defaults to False.
        show_imbalance (bool, optional): If unbalanced reaction is found, show
            output. Defaults to True.
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathway are not available.
            Defaults to: "universal"
    """
    if not isinstance(model, Model):
        raise TypeError("Model is invalid")
    if not filename:
        filename = Path.cwd().joinpath("summary.txt")
    # Retrieve information for summary methods
    basic_info = get_basic_info(model=model)
    old_values = get_DataList(model=model)
    try:
        # From identifier
        data_dict = get_data(
            directory=directory,
            identifier=str(pathway),
            database=database,
            model_id=model_id,
        )
        _from_data(
            model=model,
            data_dict=data_dict,
            directory=directory,
            database=database,
            compartment=compartment,
            avoid_list=avoid_list,
            replacement=replacement,
            ignore_list=ignore_list,
            minimum=minimum,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            model_id=model_id,
        )
    # TODO: replace Warning for specific error
    except (HTTPError, Warning):
        # From Reaction
        _from_sequence(
            model=model,
            identifier=group,
            sequence=list(pathway),
            compartment=compartment,
            avoid_list=avoid_list,
            directory=directory,
            database=database,
            replacement=replacement,
            ignore_list=ignore_list,
            minimum=minimum,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            model_id=model_id,
        )
    except TypeError:
        raise Warning("Argument 'pathway' must be iterable or a identifier")
    finally:
        # Print summary
        check_to_write(
            model=model,
            summary=summary,
            filename=filename,
            basic_info=basic_info,
            old_values=old_values,
        )
