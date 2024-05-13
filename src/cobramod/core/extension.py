"""Pathway extension

This module handles the addition of reactions as a pathway into
a model and the corresponding test that comes with it.

Most important functions:
- add_pathway: Adds a pathway or multiple reactions into a model.
- test_non_zero_flux: Checks that the given reaction in a model is active and
    gives a non-zero flux.
"""

from __future__ import annotations

from contextlib import suppress
from pathlib import Path
from typing import Generator, Optional, Union

import cobra.core as cobra_core
import cobra.exceptions as cobra_exceptions

import cobramod.retrieval as cmod_retrieval
import cobramod.utils as cmod_utils
from cobramod.core import creation as cmod_core_creation
from cobramod.core import graph as cmod_core_graph
from cobramod.core.pathway import Pathway
from cobramod.core.summary import DataModel
from cobramod.core.summary import summary as summarize
from cobramod.debug import debug_log

# Defaults to 1E-07
cobra_tolerance: float = cobra_core.Configuration().tolerance


def yield_reaction_from_list(
    sequence: list[str],
    compartment: str,
    directory: Path,
    database: Optional[str],
    show_imbalance: bool,
    stop_imbalance: bool,
    replacement: dict,
    model: cobra_core.Model,
    model_id: str,
    genome: Optional[str],
) -> Generator[cobra_core.Reaction, None, None]:
    """
    Yields a Reaction from given list of identifiers .

    .. hint:: Hyphens will become underscores. Double hyphens become single\
    underscores.

    Args:
        sequence (list): Identifiers for the reactions.
        directory (Path): Path to directory where data is located.
        database (Optional[str]): Name of the database. Check
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
        try:
            obj = model.reactions.get_by_id(identifier)
        except KeyError:
            # Check if translation is available
            obj = cmod_core_creation.create_object(
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
        if not isinstance(obj, cobra_core.Reaction):
            raise TypeError("Given object is not a valid COBRApy Reaction")
        yield obj


def find_problem(model: cobra_core.Model, identifier: str) -> list[str]:
    """
    Return a List with the metabolite identifiers that must have a sink
    reaction to make the reaction carry a non-zero flux.
    """
    problem: list[str] = list()
    value: Optional[float]

    reaction = model.reactions.get_by_id(identifier)
    if not isinstance(reaction, cobra_core.Reaction):
        raise TypeError("Given object is not a COBRApy Reaction")

    metabolites: list[cobra_core.Metabolite] = (
        reaction.products + reaction.reactants
    )

    # Add sink
    for metabolite in metabolites:
        with suppress(ValueError):
            model.add_boundary(metabolite=metabolite, type="sink")
            debug_log.debug(f"Sink reaction for {identifier} created.")

    assert model.slim_optimize() != 0.0

    # Test each reaction
    for metabolite in metabolites:
        sink = model.reactions.get_by_id(f"SK_{metabolite.id}")
        if not isinstance(sink, cobra_core.Reaction):
            raise TypeError("Given object is not a COBRApy Reaction")

        if len(metabolite.reactions) < 2:
            # Sink cannot be removed because it needs 2 reaction
            continue

        model.remove_reactions([sink])
        debug_log.debug(f"Sink reaction for {identifier} removed.")

        # Find metabolites that are necessary to carry a flux
        try:
            value = model.slim_optimize(error_value=None)

            if not value:
                value = 0.0

            if abs(value) <= cobra_tolerance:
                raise cobra_exceptions.OptimizationError("FAILED")

        except cobra_exceptions.OptimizationError:
            model.add_reactions(reaction_list=[sink])
            debug_log.debug(
                f"Sink reaction for {identifier} added again. Metabolite "
                f'"{metabolite.id}" must have a sink reaction to make '
                f'"{identifier}" carry a flux.'
            )
            assert metabolite.id
            problem.append(metabolite.id)

    return problem


def recursive_flux_test(
    model: cobra_core.Model, identifier: str, times: int = 0
) -> float:
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
        debug_log.debug(
            f'Test to carry non-zero flux for "{identifier}" started.'
        )
    # run
    passed: bool
    value = model.slim_optimize()

    if not value:
        passed = False
        value = 0.0

    elif abs(value) > cobra_tolerance:
        passed = True

    else:
        passed = False

    # Add sinks and get new value
    if not passed:
        if times < 1:
            cmod_utils.confirm_sink(model=model, identifier=identifier)
            value = recursive_flux_test(
                model=model, identifier=identifier, times=times + 1
            )

        else:
            problem = find_problem(model=model, identifier=identifier)
            msg = (
                "The model cannot turnover the following metabolites "
                f"{problem}. To overcome this, sink reactions were created to "
                "simulate their synthesis."
            )
            debug_log.warning(msg)
            value = model.slim_optimize()

    # Raise only if manual intervention is necessary
    if times == 0:
        if abs(value) > cobra_tolerance:
            debug_log.info(
                f"Non-zero flux test for reaction '{identifier}' passed."
            )

        else:
            raise cobra_exceptions.OptimizationError(
                f"There is a problem with reaction {identifier}. Please create "
                "manually the corresponding turnover reactions for the "
                f'metabolites of "{identifier}".'
            )
    return value


def non_zero_core(model: cobra_core.Model, identifier: str):
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

    reaction = model.reactions.get_by_id(identifier)

    # check reversibility
    if cmod_utils.reaction_is_minimize(model=model, identifier=identifier):
        model.objective_direction = "min"

    # Define new objective
    model.objective = reaction

    # Recursive function
    recursive_flux_test(model=model, identifier=identifier)

    # Revert back objective
    model.objective = original_objective
    model.objective_direction = original_direction


def test_non_zero_flux(model: cobra_core.Model, reaction: str):
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
    previous_sinks: set[str] = {sink.id for sink in model.sinks if sink.id}

    non_zero_core(model=model, identifier=reaction)

    # Revert back objective and inform
    cmod_utils.inform_new_sinks(model=model, previous_sinks=previous_sinks)
    print(f"Reaction {reaction} passed the non-zero-flux test.")


def add_reactions_to_Pathway(
    model: cobra_core.Model,
    pathway: Union[cobra_core.Group, Pathway],
    sequence: list[cobra_core.Reaction],
    ignore_list: list,
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
    if not all(
        (isinstance(reaction, cobra_core.Reaction) for reaction in sequence)
    ):
        raise TypeError("Reactions are not valid objects. Check list.")

    # Add sequence to model
    for reaction in sequence:
        # Add reaction if not in model
        if reaction not in model.reactions:
            model.add_reactions([reaction])
            debug_log.info(f'Reaction "{reaction.id}" was added to model.')

        # Skip test but include reaction in pathway
        if reaction.id not in ignore_list:
            non_zero_core(model=model, identifier=reaction.id)

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
    if not model.groups.has_id(pathway.id) and len(pathway.members) > 0:
        model.add_groups(group_list=[pathway])
        debug_log.info(f'Pathway "{pathway.id}" added to Model.')


def remove_avoid_reactions(
    sequence: list[str], avoid_list: list[str]
) -> list[str]:
    """
    Updates the given sequence taking into consideration reactions to avoid
    in the sequence

    Args:
        sequence (list): identifier of the reactions to update
        avoid_list (list): Reactions that should not be included in the model
    """
    new_sequence: list[str] = list()

    # Remove reactions
    for item in sequence:
        if item in avoid_list:
            msg = (
                f'Reaction "{item}" was found in "avoid_list". Reaction {item}'
                + " will not be added to the model."
            )
            debug_log.warning(msg=msg)
            continue

        new_sequence.append(item)

    return new_sequence


def add_pathway_from_file(
    model: cobra_core.Model,
    file: Union[str, Path],
    identifier: str,
    database: Optional[str],
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

     - Path or str: A file with components. E. g. :
        Path.cwd().joinpath("file_with_names.txt") or "./file_with_names.txt"

      This applies for both options :

        :code:`reaction_identifier, compartment`

        For custom reactions

        :code:`reaction_identifier, reaction_name | coefficient metabolite <->
        coefficient metabolite

    Args:
        model (Model): Model to expand.
        path (Union[str, Path]): Location of the file.
        database (str): Name of the database.
            Check :obj:`cobramod.available_databases` for a list of names.
        identifier (str): Common :class:`cobramod.pathway.Pathway`
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
    if isinstance(file, str):
        file = Path(file).absolute()

    # Either create a Pathway or give new name.
    try:
        pathway = model.groups.get_by_id(identifier)
        if not isinstance(pathway, cobra_core.Group) or not isinstance(
            pathway, Pathway
        ):
            raise TypeError(
                "Given object is not a valid COBRApy Group or a Cobramod Pathway"
            )
    except KeyError:
        pathway = Pathway(id=identifier)
        model.add_groups([pathway])

    # Get sinks
    previous_sinks: set[str] = {sink.id for sink in model.sinks if sink.id}

    # Get reactions from file
    new_reactions = cmod_core_creation.get_file_reactions(
        model=model,
        filename=file,
        directory=directory,
        database=database,
        replacement=replacement,
        stop_imbalance=stop_imbalance,
        show_imbalance=show_imbalance,
        genome=genome,
        model_id=model_id,
    )

    # Inform about sinks
    cmod_utils.inform_new_sinks(model=model, previous_sinks=previous_sinks)

    graph = cmod_core_graph.build_lineal_graph(
        sequence=[reaction.id for reaction in new_reactions]
    )

    add_reactions_to_Pathway(
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


def add_pathway_from_data(
    model: cobra_core.Model,
    data: cmod_retrieval.Data,
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
    mapping = cmod_core_graph.get_graph_dict(data.attributes["pathway"])

    # Either create a Pathway or obtain the correct Pathway.
    identifier = data.identifier

    if group:
        identifier = group

    try:
        pathway = model.groups.get_by_id(identifier)
        if not isinstance(pathway, cobra_core.Group) or not isinstance(
            pathway, Pathway
        ):
            raise TypeError(
                "Given object is not a valid COBRApy Group or a Cobramod Pathway"
            )

    except KeyError:
        pathway = Pathway(id=identifier)
        model.add_groups([pathway])

    # Get sinks
    previous_sinks: set[str] = {sink.id for sink in model.sinks if sink.id}

    # FIXME: this should be changed
    graph: dict[str, Union[tuple[str], str, None]] = dict()

    for sequence in mapping:
        # Update sequence with removal of reactions
        sequence = remove_avoid_reactions(
            sequence=sequence, avoid_list=avoid_list
        )

        if not sequence:
            continue

        # Replacements will be used here
        sequence = list(
            yield_reaction_from_list(
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
        previous: str = ""
        for i, reaction in enumerate(sequence):
            if i == len(sequence) - 1:
                graph[reaction.id] = None  # type: ignore

            if previous:
                value = graph.get(previous)

                if not value:
                    value = reaction.id  # type: ignore

                elif isinstance(value, tuple):
                    a = list(value) + [reaction.id]  # type: ignore
                    value = tuple(a)  # type: ignore

                graph[previous] = value

            previous = reaction.id  # type: ignore

        # Add to model
        add_reactions_to_Pathway(
            model=model,
            pathway=pathway,
            sequence=sequence,
            ignore_list=ignore_list,
        )

    # Inform about sinks
    cmod_utils.inform_new_sinks(model=model, previous_sinks=previous_sinks)

    if not pathway.graph:
        pathway.graph = graph

    else:
        pathway.graph.update(graph)
    pathway.notes["ORDER"] = pathway.graph


def add_pathway_from_strings(
    model: cobra_core.Model,
    identifier: str,
    sequence: list[str],
    database: Optional[str],
    compartment: str,
    directory: Path,
    avoid_list: list[str],
    replacement: dict[str, str],
    ignore_list: list[str],
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
        database (Optional[str]): Name of the database.
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
        if not isinstance(pathway, cobra_core.Group):
            raise TypeError(
                f"The group for {identifier} is not a COBRApy Group object"
            )

    except KeyError:
        pathway = Pathway(id=identifier)
        model.add_groups([pathway])

    # Data storage is managed by the function
    sequence = remove_avoid_reactions(sequence=sequence, avoid_list=avoid_list)

    # Get sinks
    previous_sinks: set[str] = {sink.id for sink in model.sinks if sink.id}

    reactions = list(
        yield_reaction_from_list(
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
    # Create graph. It will be always simple lineal
    graph = cmod_core_graph.build_lineal_graph(
        sequence=[reaction.id for reaction in reactions]
    )

    # Append to model
    add_reactions_to_Pathway(
        model=model,
        sequence=reactions,
        pathway=pathway,
        ignore_list=ignore_list,
    )

    # Inform about sinks
    cmod_utils.inform_new_sinks(model=model, previous_sinks=previous_sinks)

    if isinstance(pathway, Pathway):
        if not pathway.graph:
            pathway.graph = graph
        else:
            pathway.graph.update(graph)
        pathway.notes["ORDER"] = pathway.graph


def add_pathway(
    model: cobra_core.Model,
    pathway: Union[list[str], str, Path],
    directory: Union[Path, str],
    compartment: str,
    database: Optional[str] = None,
    group: Optional[str] = None,
    avoid_list: list[str] = [],
    replacement: dict = {},
    ignore_list: list[str] = [],
    filename: Optional[Union[str, Path]] = None,
    stop_imbalance: bool = False,
    show_imbalance: bool = True,
    model_id: str = "",
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
    if not isinstance(model, cobra_core.Model):
        raise TypeError("Model is invalid")

    if isinstance(directory, str):
        directory = Path(directory).absolute()

    if not directory.exists():
        raise FileNotFoundError(
            f"Directory '{str(directory)}' not found. Create the data directory"
        )

    # Save information for summary methods
    old_values = DataModel.from_model(model)

    # Check if identifier
    if isinstance(pathway, str):
        identifier = pathway

        pathway = Path(pathway).absolute()

        # Identifier found
        if not pathway.suffix:
            pathway = identifier

    if isinstance(pathway, str):
        data_dict = cmod_retrieval.get_data(
            directory=directory,
            identifier=pathway,
            database=database,
            model_id=model_id,
            genome=genome,
        )
        if not database:
            raise AttributeError(
                "Database argument cannot be empty. Specify the name"
            )
        add_pathway_from_data(
            model=model,
            data=data_dict,
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

    elif isinstance(pathway, Path):
        if not group:
            group = "custom_group"

        add_pathway_from_file(
            model=model,
            file=pathway,
            database=database,
            ignore_list=ignore_list,
            genome=genome,
            model_id=model_id,
            directory=directory,
            replacement=replacement,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            identifier=group,
        )

    elif isinstance(pathway, list):
        if not group:
            group = "custom_group"

        # if not database:
        #     raise AttributeError(
        #         "Database argument cannot be empty. Specify the name"
        #     )

        add_pathway_from_strings(
            model=model,
            identifier=group,
            sequence=pathway,
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

    summarize(model, old_values, filename=filename)
