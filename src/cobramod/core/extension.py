#!/usr/bin/env python3
"""Pathway extension

This module handles the addition of reactions as a pathway into
a model and the corresponding test that comes with it.

Most important functions:
- add_pathway: Adds a pathway or multiple reactions into a model.
- test_non_zero_flux: Checks that the given reaction in a model is active and
    gives a non-zero flux.
"""
from contextlib import suppress
from pathlib import Path
from typing import Union, Generator, Optional, List
from warnings import warn

from cobra import Model, Reaction
from cobra.core.configuration import Configuration
from cobra.exceptions import OptimizationError

from cobramod.core.creation import create_object
from cobramod.core.graph import build_graph, _create_quick_graph, _format_graph
from cobramod.core.pathway import Pathway
from cobramod.core.retrieval import get_data
from cobramod.core.summary import DataModel, summary as summarize
from cobramod.debug import debug_log
from cobramod.error import NotInRangeError


# Defaults to 1E-07
cobra_tolerance: float = Configuration().tolerance


def _create_reactions(
    sequence: List[str],
    compartment: str,
    directory: Path,
    database: str,
    show_imbalance: bool,
    stop_imbalance: bool,
    replacement: dict,
    model: Model,
    model_id: str,
    genome: Optional[str],
) -> Generator:
    """
    Yields a Reaction from given list of identifiers .

    .. hint:: Hyphens will become underscores. Double hyphens become single\
    underscores.

    Args:
        sequence (list): Identifiers for the reactions.
        directory (Path): Path to directory where data is located.
        database (str): Name of the database. Check
            :obj:`cobramod.available_databases` for a list of names.
        compartment (str): Location of the reaction.
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers. This applies to metabolites.
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        model_id (str): Exclusive for BIGG. Retrieve object from specified
            model.
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained for this specie.

    Returns:
        Generator: New Reactions objects
    """

    debug_log.debug(f"Obtaining data for following reactions {sequence}.")
    # From given list (which could include None values), retrieve only Reaction
    # Objects.
    for identifier in sequence:
        # Check if translation is available
        obj = create_object(
            identifier=identifier,
            directory=directory,
            database=database,
            compartment=compartment,
            replacement=replacement,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            model=model,
            model_id=model_id,
            genome=genome,
        )
        yield obj


def _find_next_demand(
    model: Model, reaction_id: str, ignore_list: List[str] = []
) -> str:
    """
    Returns an identifier of first metabolite found either in the product or
    reactant side of given reaction.

    Reversibility of the reaction is taken into consideration. A list with
    metabolites identifiers can be passed to ignored them.

    Args:
        model (Model): model to test
        reaction_id (str): reaction identifier for model
        ignore_list (list, optional): A sequence of formatted metabolites
            to ignore. Defaults to []

    Raises:
        Warning: if no metabolite was found

    Returns:
        str: metabolite identifier to create a demand reaction
    """
    reaction: Reaction = model.reactions.get_by_id(reaction_id)
    lb, up = reaction.bounds
    # Checking reversibility
    metabolites = {
        "products": (
            metabolite.id
            for metabolite in reaction.products
            if metabolite.id not in ignore_list
        ),
        "reactants": (
            metabolite.id
            for metabolite in reaction.reactants
            if metabolite.id not in ignore_list
        ),
    }
    # left --> right
    if up > abs(lb):
        candidates = metabolites["products"]
    # left <-- right
    elif up < abs(lb):
        candidates = metabolites["reactants"]
    # left <--> right
    else:
        # TODO: decide what to do (reversibility)
        # FIXME: isomers sometimes shows double demand
        candidates = metabolites["products"]
    try:
        demand: str = next(candidates)
        return demand
    except StopIteration:
        raise StopIteration(
            "No metabolite found to become a demand for the non-zero flux test"
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
    reaction = f"{type_prefix[boundary]}{metabolite}"
    if reaction in (
        reaction.id
        for reaction in model.metabolites.get_by_id(metabolite).reactions
    ):
        model.remove_reactions([reaction])
        debug_log.warning(
            f'Auxiliary {boundary} reaction for "{metabolite}" removed.'
        )


def _verify_boundary(model: Model, metabolite: str, ignore_list: List[str]):
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
        ignore_list (list): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []

    Raises:
        Warning: If a metabolite is found in given ignore list.
    """
    debug_log.debug(
        "Counting number of auxiliary demand and sink reactions related to "
        f'metabolite "{metabolite}".'
    )
    if metabolite in ignore_list:
        msg = (
            f'Metabolite "{metabolite}" was found in "ignore_list". No '
            "auxiliary sink reactions will be be created during the non-zero "
            " flux test."
        )
        debug_log.warning(msg)
        raise Warning(msg)
    reactions = len(model.metabolites.get_by_id(metabolite).reactions)
    demands = [
        reaction
        for reaction in model.metabolites.get_by_id(metabolite).reactions
        if "DM_" in reaction.id
    ]
    # Check to add sinks. If demand reaction is found and metabolite has
    # less than two reactions, then create sink, otherwise remove any sink.
    if demands:
        if reactions <= 2:
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite), type="sink"
            )
            debug_log.debug(
                msg=f'Auxiliary sink reaction for "{metabolite}" created.'
            )

        else:
            _remove_boundary(
                model=model, metabolite=metabolite, boundary="sink"
            )
    # Otherwise, if it has no demand, add the sink if it has only one reaction
    # and remove sink if more than 2 reactions.
    else:
        if reactions <= 1:
            model.add_boundary(
                metabolite=model.metabolites.get_by_id(metabolite), type="sink"
            )
            debug_log.debug(
                msg=f'Auxiliary sink reaction for "{metabolite}" created.'
            )
        elif not reactions <= 2:
            _remove_boundary(
                model=model, metabolite=metabolite, boundary="sink"
            )
    # Remove demand if needed. This part is always necessary since the extra
    # demands would add extra flux into the system.
    _remove_boundary(model=model, metabolite=metabolite, boundary="demand")


def _fix_side(
    model: Model, reaction: str, side: str, ignore_list: List[str] = []
):
    """
    Checks for either the product or reactant side of a reactions, if
    participant-metabolites have enough sink reactions to produce a feasible
    answer. If necessary, it will creates or remove them.

    Args:
        model (Model): model to test.
        reaction_id (str): reaction identifier for given model.
        side (str): Side to verify. Options: "right", "left"
        ignore_list (list, optional): A sequence of formatted metabolites
            to ignore when testing new reactions. Defaults to []:

    Raises:
        ValueError: if side is not in the options
    """
    participants = {
        "right": model.reactions.get_by_id(reaction).products,
        "left": model.reactions.get_by_id(reaction).reactants,
    }
    try:
        metabolites = participants[side]
    except KeyError:
        raise ValueError('Only valid options are "right" and "left"')
    for meta in metabolites:
        # create or remove
        with suppress(Warning):
            _verify_boundary(
                model=model, metabolite=meta.id, ignore_list=ignore_list
            )


def _verify_sinks(model: Model, reaction: str, ignore_list: List[str]):
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
        ignore_list (list): A sequence of formatted metabolites
            to ignore when testing new reactions.
    """
    debug_log.debug(
        "Counting number of auxiliary demand and sink reactions for reaction "
        f'"{reaction}". Auxiliary reactions ensure that all reactants and '
        "products of the equation can be turned over in the model."
    )
    # reactant side
    _fix_side(
        model=model, reaction=reaction, ignore_list=ignore_list, side="left"
    )
    # product side
    _fix_side(
        model=model, reaction=reaction, ignore_list=ignore_list, side="right"
    )


def test_non_zero_flux(
    model: Model, reaction: str, times: int = 0, ignore_list: List[str] = []
):
    """
    Verifies that a simple FBA can be performed. A demand reaction is created
    with a minimum flux based on the COBRApy Configuration object. It will use
    its variable 'tolerance' multiplied by 10. The function is recursive and
    checks if sink reactions are sufficient or exceeded. It generates a demand
    reaction for the reaction used for testing and removes it if necessary.

    Args:
        model (Model): Model in which the reaction is present.
        reaction (str): Reaction identifier in model to test.
        times (int, optional): Track of recursions. Defaults to 0.
        ignore_list (list, optional): A sequence of formatted metabolites to be
        ignored when testing new reactions. Defaults to []

    Raises:
        NotInRangeError: If the solution is infeasible after many recursions.
            Depends on the number of metabolites in the reaction.

    """
    if times == 0:
        debug_log.info(
            f'Test to carry non-zero fluxes for "{reaction}" started'
        )
    # finding demand for testing
    next_demand = _find_next_demand(
        model=model, reaction_id=reaction, ignore_list=ignore_list
    )
    # The first step is to add a demand for check that the reaction has an
    # active flux.
    with suppress(ValueError):
        model.add_boundary(model.metabolites.get_by_id(next_demand), "demand")
        model.reactions.get_by_id(f"DM_{next_demand}").lower_bound = abs(
            cobra_tolerance * 10
        )
        debug_log.debug(
            f'Demand reaction "DM_{next_demand}" for non-zero flux test was '
            + "added to the model."
        )
    # Setting maximum times for recursion
    if times == len(model.reactions.get_by_id(reaction).metabolites):
        raise NotInRangeError(reaction=reaction)
    # answer must be reasonable
    # comparison must be using absolute values
    try:
        model.slim_optimize(error_value=None)
        # if works, pass and return old objective
        debug_log.info(f'Reaction "{reaction}" passed the non-zero flux test.')
        # Remove old demand
        _remove_boundary(
            model=model, metabolite=next_demand, boundary="demand"
        )
        _verify_sinks(model=model, reaction=reaction, ignore_list=ignore_list)
    except OptimizationError:
        debug_log.debug(
            f'Non-zero flux test for reaction "{reaction}" failed. Flux value'
            f" lower than solver tolerance ({cobra_tolerance*10}) defined in "
            " cobramod.core.extension.solver_tolerance}. Please check the"
            " logging information."
        )
        _verify_sinks(model=model, reaction=reaction, ignore_list=ignore_list)
        # Recursive with 'extra' argument
        test_non_zero_flux(
            model=model,
            reaction=reaction,
            times=times + 1,
            ignore_list=ignore_list,
        )


def _add_sequence(
    model: Model, pathway: Pathway, sequence: list, ignore_list: list
):
    """
    From a sequence of Reaction objects, add each Reaction into given model. It
    checks if new reactions do not break the optimized value. All reactions
    will be added to a common COBRApy Group.

    Args:
        model (Model): Model to expand.
        pathway (Pathway): Common Pathway to add the reaction.
        sequence (list): List with :class:`cobra.core.reaction.Reaction`
            objects
        ignore_list (list): A sequence of formatted metabolites to skip when
            testing, and/or reactions that should be added but not tested.
            This is useful for long cyclical pathways.

    Raises:
        TypeError: if reactions are not valid Reaction objects
    """
    if not all((isinstance(reaction, Reaction) for reaction in sequence)):
        raise TypeError("Reactions are not valid objects. Check list.")
    # Add sequence to model
    reaction: Reaction
    for reaction in sequence:
        # Add reaction if not in model
        if reaction not in model.reactions:
            model.add_reactions([reaction])
            debug_log.info(f'Reaction "{reaction.id}" was added to model.')
        # Skip test but include reaction in pathway
        if reaction.id not in ignore_list:
            test_non_zero_flux(
                model=model, reaction=reaction.id, ignore_list=ignore_list
            )
        else:
            debug_log.warning(
                f'Reaction "{reaction.id}" found in "ignore_list". Skipping '
                + "non-zero flux test."
            )
            # Add to pathway only if reaction was not previously in the model.
        debug_log.debug(
            f'Reaction "{reaction.id}" passed the non-zero flux test.'
        )
        pathway.add_members(new_members=[reaction])
        debug_log.info(
            f'Reaction "{reaction.id}" added to group "{pathway.id}".'
        )
    debug_log.debug(f'Reactions added to group "{pathway.id}".')
    # Only add if there is at least 1 reaction in the group.
    if (
        pathway.id not in [group.id for group in model.groups]
        and len(pathway.members) > 0
    ):
        model.add_groups(group_list=[pathway])
        debug_log.info(f'Pathway "{pathway.id}" added to Model.')


def _update_reactions(
    sequence: List[str], avoid_list: List[str], replacement: dict
) -> List[str]:
    """
    Updates the given sequence taking into consideration reactions to avoid
    and to rename

    Args:
        sequence (list): Names of the reactions to update
        avoid_list (list): Reactions that should not be included in the model
        replacement (dict): Directory with identifiers to be replaced. Keys are
            the original identifiers and the values the new identifiers
    """
    new_sequence: List[str] = list()
    for item in sequence:
        # Remove reactions
        if item in avoid_list:
            msg = (
                f'Reaction "{item}" was found in "avoid_list". Reaction {item}'
                + " will not be added to the model."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
            continue
        # Replace reactions
        if item in replacement.keys():
            msg = (
                f'Reaction "{item}" was found in "replacement". Reaction '
                f'will be replaced by "{replacement[item]}".'
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
            item = replacement[item]
        # add in the end
        new_sequence.append(item)
    return new_sequence


def _from_data(
    model: Model,
    data_dict: dict,
    group: Optional[str],
    directory: Path,
    database: str,
    compartment: str,
    avoid_list: list,
    replacement: dict,
    ignore_list: list,
    stop_imbalance: bool,
    show_imbalance: bool,
    model_id: str,
    genome: Optional[str],
):
    """
    Adds a pathway into given model from a dictionary with the information of
    the pathway into given model from a dictionary with the information of
    the model.

    Args:
        model (Model): Model to expand.
        data_dict (dict): Dictionary with the information for the pathway.
        directory (Path): Path for directory to stored and retrieve data.
        database (str): Name of the database to search for reactions and
            metabolites. Check :obj:`cobramod.available_databases` for a list
            of names.
        compartment: Location of the reactions.

    Arguments for complex pathways:
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
        avoid_list (list): A sequence of reactions identifiers to
            avoid adding to the model. This is useful for long pathways,
            where X reactions need to be excluded.
        ignore_list (list): A sequence of formatted metabolites to skip when
            testing, and/or reactions that should be added but not tested.
            This is useful for long cyclical pathways.

    Arguments for utilities:
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathway are not available.
            Defaults to: "universal"
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie.
            List available at https://www.genome.jp/kegg/catalog/org_list.html
    """
    # Create mapping from dictionary
    mapping = build_graph(graph=data_dict["PATHWAY"])
    # Either create a Pathway or obtain the correct Pathway.
    identifier = data_dict["ENTRY"]
    if group:
        identifier = group
    try:
        pathway = model.groups.get_by_id(identifier)
    except KeyError:
        pathway = Pathway(id=identifier)
    # Get sinks
    previous_sink = {sink.id for sink in model.sinks}
    for sequence in mapping:
        # Update sequence with new identifiers or removal of reactions
        sequence = _update_reactions(
            sequence=sequence, avoid_list=avoid_list, replacement=replacement
        )
        if not sequence:
            continue
        # Data storage is handled by method. Reactions objects builder:
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
                genome=genome,
            )
        )
        # Add to model
        _add_sequence(
            model=model,
            pathway=pathway,
            sequence=sequence,
            ignore_list=ignore_list,
        )
    # Inform about sinks
    sinks = {sink.id for sink in model.sinks}.difference(previous_sink)
    if sinks:
        for reaction in sinks:
            msg = (
                f'Auxiliary sink reaction for "{reaction}" created. Consider '
                "removing it and adding the synthesis reactions for the "
                "metabolite."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)

    # Add graph to Pathway
    graph = _format_graph(
        graph=data_dict["PATHWAY"],
        model=model,
        compartment=compartment,
        database=database,
        model_id=model_id,
        directory=directory,
        avoid_list=avoid_list,
        replacement=replacement,
        genome=genome,
    )
    if not pathway.graph:
        pathway.graph = graph
    else:
        pathway.graph.update(graph)
    pathway.notes["ORDER"] = pathway.graph


def _from_sequence(
    model: Model,
    identifier: str,
    sequence: List[str],
    database: str,
    compartment: str,
    directory: Path,
    avoid_list: List[str],
    replacement: dict,
    ignore_list: List[str],
    stop_imbalance: bool,
    show_imbalance: bool,
    model_id: str,
    genome: Optional[str],
):
    """
    Adds a sequence of identifiers to given model. It will automatically test
    for valid optimized values.

    Args:
        model (Model): Model to expand.
        identifier (str): Common :class:`cobramod.pathway.Pathway` identifier.
        sequence (list): List reaction identifiers.
        database (str): Name of the database.
            Check :obj:`cobramod.available_databases` for a list of names.
        compartment: Location of the reactions.
        directory (Path): Path for directory to stored and retrieve data.

    Arguments for complex pathways:
        avoid_list (list): A sequence of reactions identifiers to
            avoid adding to the model. This is useful for long pathways,
            where X reactions need to be excluded.
        ignore_list (list): A sequence of formatted metabolites to skip when
            testing, and/or reactions that should be added but not tested.
            This is useful for long cyclical pathways.
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.

    Arguments for utilities:
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathway are not available.
        genome (str): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie.
            List available at https://www.genome.jp/kegg/catalog/org_list.html
    """
    # Either create a Pathway or obtain the correct Pathway.
    try:
        pathway = model.groups.get_by_id(identifier)
    except KeyError:
        pathway = Pathway(id=identifier)
    # Create graph. It will be always simple lineal
    graph = _create_quick_graph(sequence=sequence)
    # Data storage is managed by the function
    sequence = _update_reactions(
        sequence=sequence, avoid_list=avoid_list, replacement=replacement
    )
    # Get sinks
    previous_sink = {sink.id for sink in model.sinks}
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
            genome=genome,
        )
    )
    # Append to model
    _add_sequence(
        model=model,
        sequence=sequence,
        pathway=pathway,
        ignore_list=ignore_list,
    )
    # Inform about sinks
    sinks = {sink.id for sink in model.sinks}.difference(previous_sink)
    if sinks:
        for reaction in sinks:
            msg = (
                f'Auxiliary sink reaction for "{reaction}" created. Consider '
                "removing it and adding the synthesis reactions for the "
                "metabolite."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
    # Add graph to Pathway
    graph = _format_graph(
        graph=graph,
        model=model,
        compartment=compartment,
        database=database,
        model_id=model_id,
        directory=directory,
        avoid_list=avoid_list,
        replacement=replacement,
        genome=genome,
    )
    # FIXME: pointers?
    if not pathway.graph:
        pathway.graph = graph
    else:
        pathway.graph.update(graph)
    pathway.notes["ORDER"] = pathway.graph


def add_pathway(
    model: Model,
    pathway: Union[list, str],
    directory: Path,
    database: str,
    compartment: str,
    group: Optional[str] = None,
    avoid_list: List[str] = [],
    replacement: dict = {},
    ignore_list: List[str] = [],
    filename: Path = None,
    stop_imbalance: bool = False,
    show_imbalance: bool = True,
    model_id: str = "universal",
    genome: Optional[str] = None,
):
    """
    Adds a pathway from the given database into a model. The argument 'pathway'
    can be a list of reactions identifiers or explicitly a pathway identifier.
    A group of reactions will be included in a custom group. The data will be
    downloaded and structured according to the database.

    Args:
        model (Model): Model to expand.
        pathway (list, str): Sequence of reaction identifiers or a single
            identifier for a pathway. Examples: ["RXN-2206", "RXN-207"] or
            "PWY-886"
        directory (Path): Path for the directory to stored and retrieve data.
        database (str): Name of the database.
            Check :obj:`cobramod.available_databases` for a list of names.
        compartment: Location of the reactions.
        group (str, optional): Common :class:`cobramod.pathway.Pathway`
            identifier. This will overwrite the name of the pathway.

    Arguments for complex pathways:
    avoid_list (list, optional): A sequence of reactions identifiers to avoid
        adding to the model. This is useful for long pathways, where
        X reactions need to be excluded.
    replacement (dict, optional): Original identifiers to be replaced.
        Values are the new identifiers.
    ignore_list (list): A sequence of formatted metabolites to be ignored when
    testing, and/or reactions that should be added but not tested.
        This is useful for long cyclical pathways.

    Arguments for summary:
        filename (Path, optional): Location for the summary. Defaults to
            "summary" in the current working directory. The file format is
            defined by the suffix. The suffixes '.txt', '.csv' and '.xlsx' can
            be used. If the filename is set to None, no summary will be
            created.

    Arguments for utilities:
        stop_imbalance (bool, optional): If an unbalanced reaction is found,
            stop the process. Defaults to False.
        show_imbalance (bool, optional): If an unbalanced reaction is found,
            print output. Defaults to True.
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathways are not available.
            Defaults to: "universal"
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            species involved. Genes will be obtained for this species.
            List available at https://www.genome.jp/kegg/catalog/org_list.html
    """
    if not isinstance(model, Model):
        raise TypeError("Model is invalid")
    # Save information for summary methods
    old_values = DataModel.from_model(model)
    if isinstance(pathway, str):
        # Get data and transform to a pathway
        data_dict = get_data(
            directory=directory,
            identifier=str(pathway),
            database=database,
            model_id=model_id,
            genome=genome,
        )
        # Run the function to convert the reaction, create the graph and add
        # to Pathway
        _from_data(
            model=model,
            data_dict=data_dict,
            directory=directory,
            database=database,
            compartment=compartment,
            avoid_list=avoid_list,
            replacement=replacement,
            ignore_list=ignore_list,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            model_id=model_id,
            genome=genome,
            group=group,
        )
    elif isinstance(pathway, list):
        # From Reaction
        if not group:
            group = "custom_group"
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
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            model_id=model_id,
            genome=genome,
        )
    else:
        raise ValueError("Argument 'pathway' must be iterable or a identifier")
    # Print summary
    summarize(model, old_values, filename=filename)
