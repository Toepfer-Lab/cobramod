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
from typing import Generator, List, Optional, Set, Union
from warnings import warn

from cobra import Metabolite, Model, Reaction
from cobra.core.configuration import Configuration
from cobra.exceptions import OptimizationError

from cobramod.core.creation import (
    _convert_string_reaction,
    _get_file_reactions,
    create_object,
)
from cobramod.core.graph import _create_quick_graph, _format_graph, build_graph
from cobramod.core.pathway import Pathway
from cobramod.core.retrieval import get_data
from cobramod.core.summary import DataModel, summary as summarize
from cobramod.debug import debug_log


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


def _find_problem(model: Model, reaction: str) -> List[str]:
    """
    Return a List with the metabolite identifiers that must have a sink
    reaction to make the reaction carry a non-zero flux.
    """
    problem: List[str] = list()
    value: Optional[float]
    _reaction: Reaction = model.reactions.get_by_id(reaction)
    metabolites: List[Metabolite] = _reaction.products + _reaction.reactants

    # Add sink
    for metabolite in metabolites:
        with suppress(ValueError):
            model.add_boundary(metabolite=metabolite, type="sink")
            debug_log.debug(f"Sink reaction for {reaction} created.")

    assert model.slim_optimize() != 0.0

    # Test each reaction
    for metabolite in metabolites:
        sink: Reaction = model.reactions.get_by_id(f"SK_{metabolite.id}")

        if len(metabolite.reactions) < 2:
            # Sink cannot be removed because it needs 2 reaction
            continue

        model.remove_reactions([sink])
        debug_log.debug(f"Sink reaction for {reaction} removed.")

        # Find metabolites that are necessary to carry a flux
        try:
            value = model.slim_optimize(error_value=None)

            if not value:
                value == 0.0

            assert value is not None
            if abs(value) <= cobra_tolerance:
                raise OptimizationError("FAILED")

        except OptimizationError:
            model.add_reactions(reaction_list=[sink])
            debug_log.debug(
                f"Sink reaction for {reaction} added again. Metabolite "
                f'"{metabolite.id}" must have a sink reaction to make '
                f'"{reaction}" carry a flux.'
            )
            assert metabolite.id
            problem.append(metabolite.id)

    return problem


def _new_verify_sink(model: Model, reaction: str):
    """
    Creates sink reactions for the metabolites that participate in given
    reaction. They can only be created if the corresponding metabolite only
    participates in one reaction. Otherwise, the metabolite can be somehow
    turnover
    """
    _reaction: Reaction = model.reactions.get_by_id(reaction)
    metabolites: List[Metabolite] = _reaction.products + _reaction.reactants

    for metabolite in metabolites:
        amount = len(metabolite.reactions)

        if amount <= 1:
            # Warning is called in _inform_sink
            model.add_boundary(metabolite=metabolite, type="sink")


def _recursive_flux_test(
    model: Model, reaction: str, times: int = 0
) -> Optional[float]:
    """
    Main part of the non-zero flux test. It is a recursive function. Firstly,
    the function checks if there a non-zero flux. If not, it will create the
    necessary sink reactions and rerun the test. This function warns the user
    if sink reactions were created.

    Args:
        model (Model): Model where the reactions are located
        reaction (str): identifier of the reaction
        times (int): Tracks how many times this function was called.

    Returns:
        Optional[float]: optimization value

    Raises:
        OptimizationError: Whenever a reaction needs manual intervention
    """

    if times == 0:
        debug_log.info(
            f'Test to carry non-zero flux for "{reaction}" started.'
        )
    # run
    passed: bool
    value: Optional[float] = model.slim_optimize(error_value=None)

    if not value:
        passed = False
        value = 0.0

    # TODO: should be higher or at least?
    elif abs(value) > cobra_tolerance:
        passed = True

    else:
        passed = False

    # Add sinks and get new value
    if not passed:
        if times < 1:
            _new_verify_sink(model=model, reaction=reaction)
            value = _recursive_flux_test(
                model=model, reaction=reaction, times=times + 1
            )

        else:
            problem = _find_problem(model=model, reaction=reaction)
            msg = (
                "The model cannot turnover the following metabolites "
                f"{problem}. To overcome this, sink reactions were created to "
                "simulate their synthesis."
            )
            debug_log.warning(msg)
            warn(msg)
            value = model.slim_optimize()

    # Raise only if manual intervention is necessary
    assert value is not None
    if times == 0:
        if abs(value) > cobra_tolerance:
            debug_log.info(
                f"Non-zero flux test for reaction {reaction} passed."
            )

        else:
            raise OptimizationError(
                f"There is a problem with reaction {reaction}. Please create "
                "manually the corresponding turnover reactions for the "
                f'metabolites of "{reaction}".'
            )
    return value


def _is_minimize(model: Model, reaction: str) -> bool:
    """
    Return whether given reaction should be minimized when optimizing
    """
    _reaction: Reaction = model.reactions.get_by_id(reaction)
    lower: float = _reaction.lower_bound  # type: ignore
    upper: float = _reaction.upper_bound  # type: ignore

    if abs(lower) > upper:
        return True
    return False


def _test_non_zero_main(model: Model, reaction: str):
    """
    Performs non-zero flux test. In this test, a reaction is tested to make
    sure it can carry a flux. If necessary, CobraMod creates sink reactions
    to simulate turnover of metabolites that participate in the reaction.
    Warnings will be shown for these cases.

    Args:
        model (Model): Model that contains the reaction
        reaction (str): identifier of the reaction

    Raises:
        OptimizationError: If given reaction needs manual curation

    """
    # Save old objective
    original_objective = model.objective
    original_direction = model.objective_direction

    _reaction: Reaction = model.reactions.get_by_id(reaction)

    # check reversibility
    if _is_minimize(model=model, reaction=reaction):
        model.objective_direction = "min"

    # Define new objective
    model.objective = _reaction

    # Recursive function
    _recursive_flux_test(model=model, reaction=reaction)

    # Revert back objective
    model.objective = original_objective
    model.objective_direction = original_direction


def test_non_zero_flux(model: Model, reaction: str):
    """
    Performs non-zero flux test. In this test, a reaction is tested to make
    sure it can carry a flux. If necessary, CobraMod creates sink reactions
    to simulate turnover of metabolites that participate in the reaction.
    Warnings will be shown for these cases.

    Args:
        model (Model): Model that contains the reaction
        reaction (str): identifier of the reaction

    Raises:
        OptimizationError: If given reaction needs manual curation

    """
    previous_sinks: Set[str] = {sink.id for sink in model.sinks if sink.id}

    _test_non_zero_main(model=model, reaction=reaction)

    # Revert back objective and inform
    _inform_sink(model=model, previous_sinks=previous_sinks)
    print(f"Reaction {reaction} passed the non-zero-flux test.")


def _add_sequence(
    model: Model, pathway: Pathway, sequence: list, ignore_list: list
):
    """
    From a sequence of Reaction objects, add each Reaction into given model. It
    checks if new reactions do not break the optimized value. All reactions
    are added to a common CobraMod Pathway.

    Args:
        model (Model): Model to expand.
        pathway (Pathway): Common Pathway to add the reaction.
        sequence (list): List with :class:`cobra.core.reaction.Reaction`
            objects
        ignore_list (list, optional): A sequence of reactions that should be
            added but not tested for a non-zero-flux.

    Raises:
        TypeError: if reactions are not valid Reaction objects
    """
    if not all((isinstance(reaction, Reaction) for reaction in sequence)):
        raise TypeError("Reactions are not valid objects. Check list.")

    # Add sequence to model
    reaction: Reaction
    for reaction in sequence:
        assert reaction.id

        # Add reaction if not in model
        if reaction not in model.reactions:
            model.add_reactions([reaction])
            debug_log.info(f'Reaction "{reaction.id}" was added to model.')

        # Skip test but include reaction in pathway
        if reaction.id not in ignore_list:
            _test_non_zero_main(model=model, reaction=reaction.id)

        else:
            debug_log.warning(
                f'Reaction "{reaction.id}" found in "ignore_list". Skipping '
                + "non-zero flux test."
            )

        # Add to pathway only if reaction was not previously in the model.
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


def _update_reactions(sequence: List[str], avoid_list: List[str]) -> List[str]:
    """
    Updates the given sequence taking into consideration reactions to avoid
    in the sequence

    Args:
        sequence (list): identifier of the reactions to update
        avoid_list (list): Reactions that should not be included in the model
    """
    new_sequence: List[str] = list()

    # Remove reactions
    for item in sequence:
        if item in avoid_list:
            msg = (
                f'Reaction "{item}" was found in "avoid_list". Reaction {item}'
                + " will not be added to the model."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
            continue

        new_sequence.append(item)

    return new_sequence


def _inform_sink(model: Model, previous_sinks: Set[str]):
    """
    Warns the user if sinks were created during the extension of the Model.
    i.e. Function compares a set with sink identifiers with the actual sinks
    in a model.
    """
    sinks: Set[str] = {sink.id for sink in model.sinks if sink.id}.difference(
        previous_sinks
    )

    if sinks:
        for reaction in sinks:
            msg = (
                f'Auxiliary sink reaction for "{reaction}" created. Consider '
                "removing it and adding the synthesis reactions for the "
                "metabolite."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)


def _from_strings(
    model: Model,
    obj: Union[List[str], Path],
    pathway_name: str,
    database: str,
    replacement: dict,
    ignore_list: list,
    stop_imbalance: bool,
    show_imbalance: bool,
    directory: Path,
    genome: Optional[str],
    model_id: str,
):
    """
    Adds a Pathway into given model. The reactions are created from::

     - Path: A file with components. E. g. :
        Path.cwd().joinpath("file_with_names.txt")
     - List[str]: A list with multiple string. Either the identifier with its
     corresponding compartment; a string with all components or the reaction
     identifier already in the model. This applies for the Path option. E.g. :

        :code:`reaction_identifier, compartment`

        For custom reactions

        :code:`reaction_identifier, reaction_name | coefficient metabolite <->
        coefficient metabolite

    Args:
        model (Model): Model to expand.
        obj (Union[List[str], Path]): Location of the file or list with strings
        database (str): Name of the database.
            Check :obj:`cobramod.available_databases` for a list of names.
        pathway_name (str): Common :class:`cobramod.pathway.Pathway`
            identifier

    Argument for complex pathways:
        avoid_list (list, optional): A sequence of reactions identifiers to
            avoid adding to the model.
        replacement (dict, optional): Original identifiers to be replaced.
            Values are the new identifiers. This applies to metabolites as
            well. User can rename or replace identifiers using this argument.
        ignore_list (list, optional): A sequence of reactions that should be
            added but not tested for a non-zero-flux.

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
    # Either create a Pathway or give new name.
    try:
        pathway = model.groups.get_by_id(pathway_name)
    except KeyError:
        pathway = Pathway(id=pathway_name)

    # Get sinks
    previous_sinks: Set[str] = {sink.id for sink in model.sinks if sink.id}

    # Get reactions from file
    if isinstance(obj, Path):
        new_reactions = _get_file_reactions(
            model=model,
            filename=obj,
            directory=directory,
            database=database,
            replacement=replacement,
            stop_imbalance=stop_imbalance,
            show_imbalance=show_imbalance,
            genome=genome,
            model_id=model_id,
        )

    # Or convert list
    else:
        new_reactions = [
            _convert_string_reaction(
                line=line,
                model=model,
                directory=directory,
                database=database,
                replacement=replacement,
                stop_imbalance=stop_imbalance,
                show_imbalance=show_imbalance,
                genome=genome,
                model_id=model_id,
            )
            for line in obj
        ]
    # Inform about sinks
    _inform_sink(model=model, previous_sinks=previous_sinks)

    graph = _create_quick_graph(
        sequence=[reaction.id for reaction in new_reactions]
    )

    _add_sequence(
        model=model,
        pathway=pathway,
        sequence=new_reactions,
        ignore_list=ignore_list,
    )

    # No need to format graph because there shouldn't be replacements
    if not pathway.graph:
        pathway.graph = graph
    else:
        pathway.graph.update(graph)
    pathway.notes["ORDER"] = pathway.graph


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
        avoid_list (list, optional): A sequence of reactions identifiers to
            avoid adding to the model.
        replacement (dict, optional): Original identifiers to be replaced.
            Values are the new identifiers. This applies to metabolites as
            well. User can rename or replace identifiers using this argument.
        ignore_list (list, optional): A sequence of reactions that should be
            added but not tested for a non-zero-flux.

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
    previous_sinks: Set[str] = {sink.id for sink in model.sinks if sink.id}

    for sequence in mapping:
        # Update sequence with removal of reactions
        sequence = _update_reactions(sequence=sequence, avoid_list=avoid_list)
        if not sequence:
            continue

        # Replacements will be used here
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
    _inform_sink(model=model, previous_sinks=previous_sinks)

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
        directory (Path): Path for directory to store and retrieve data.

    Arguments for complex pathways:
        avoid_list (list, optional): A sequence of reactions identifiers to
            avoid adding to the model.
        replacement (dict, optional): Original identifiers to be replaced.
            Values are the new identifiers. This applies to metabolites as
            well. User can rename or replace identifiers using this argument.
        ignore_list (list, optional): A sequence of reactions that should be
            added but not tested for a non-zero-flux.

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
    sequence = _update_reactions(sequence=sequence, avoid_list=avoid_list)

    # Get sinks
    previous_sinks: Set[str] = {sink.id for sink in model.sinks if sink.id}

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
    _inform_sink(model=model, previous_sinks=previous_sinks)

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
    compartment: str,
    database: Optional[str] = None,
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
            identifier for a pathway. Examples:
            ["RXN-2206", "RXN-207"], (Reaction to download from a database)
            "PWY-886", (Pathway to download from a database)
            ["ACALDt, RXN_2206_c"], (which are reactions already in the model)
        directory (Path): Path for the directory to store and retrieve data.
        database (str, Optional): Name of the database.
            Check :obj:`cobramod.available_databases` for a list of names.
            When adding reactions of the model, this argument is not necessary
            Defaults to None.
        compartment: Location of the reactions. If adding reaction already
            in the model, this argument will not change the reaction's
            compartment.
        group (str, optional): Common :class:`cobramod.pathway.Pathway`
            identifier. This will overwrite the name of the pathway.

    Arguments for complex pathways:
        avoid_list (list, optional): A sequence of reactions identifiers to
            avoid adding to the model.
        replacement (dict, optional): Original identifiers to be replaced.
            Values are the new identifiers. This applies to metabolites as
            well. User can rename or replace identifiers using this argument.
        ignore_list (list, optional): A sequence of reactions that should be
            added but not tested for a non-zero-flux.

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

    # FIXME: static typing for database
    if isinstance(pathway, str):
        # Get data and transform to a pathway
        data_dict = get_data(
            directory=directory,
            identifier=str(pathway),
            database=database,  # type: ignore
            model_id=model_id,
            genome=genome,
        )
        # Run the function to convert the reaction, create the graph and add
        # to Pathway
        _from_data(
            model=model,
            data_dict=data_dict,
            directory=directory,
            database=database,  # type: ignore
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

        try:
            # No need for compartment and avoid_list
            _from_strings(
                model=model,
                obj=pathway,
                database=database,  # type: ignore
                ignore_list=ignore_list,
                genome=genome,
                model_id=model_id,
                directory=directory,
                replacement=replacement,
                show_imbalance=show_imbalance,
                stop_imbalance=stop_imbalance,
                pathway_name=group,
            )

        except KeyError:
            _from_sequence(
                model=model,
                identifier=group,
                sequence=list(pathway),
                compartment=compartment,
                avoid_list=avoid_list,
                directory=directory,
                database=database,  # type: ignore
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
