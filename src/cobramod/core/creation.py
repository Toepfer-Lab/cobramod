#!/usr/bin/env python3
"""Object creation

This module handles the creation of COBRApy's objects
:func:`cobra.core.metabolite.Metabolite` and
:func:`cobra.core.reaction.Reaction`. Dictionaries with the data from given
database is used. Important functions are:

- create_object: Creates and Returns a COBRApy object.
- add_reactions: Add reactions from multiple sources.
- add_metabolites: Add reactions from multiple sources

These functions are a mix of multiple simpler functions:
- _metabolite_from_string, _reaction_from_string: create objects from strings.
- _get_metabolite, _get_reaction: create objects from dictionary.
- _convert_string_reaction, _convert_string_metabolite: create objects from
files.
"""
from collections import Counter
from contextlib import suppress
from pathlib import Path
from typing import Union, Generator, List, Any

from cobra import Metabolite, Model, Reaction

from cobramod.debug import debug_log
from cobramod.error import WrongDataError, NoIntersectFound, WrongSyntax
from cobramod.core.retrieval import get_data
from cobramod.utils import _read_lines, check_imbalance, _first_item


def _build_metabolite(
    identifier: str, formula: str, name: str, charge: float, compartment: str
):
    """
    Returns a basic :funf:`cobra.Metabolite`. It will log a with a DEBUG level.

    Args:
        identifier (str): Short name for Metabolite.
        formula (str): Chemical formula
        name (name): Long name for Metabolite
        charge (float): Charge of the metabolite
        compartment (str): Location

    Returns:
        Metabolite: new created metabolite
    """
    metabolite = Metabolite(
        id=identifier,
        formula=formula,
        name=name,
        charge=charge,
        compartment=compartment,
    )
    debug_log.debug(f'Metabolite "{identifier}" created.')
    return metabolite


def _metabolite_from_string(line_string: str) -> Metabolite:
    """
    Creates a Metabolite object based on a string.

    The string must follow the syntax:

    :code:`formatted identifier, name , compartment, chemical_formula,
    molecular_charge`

    Args:
        line_string (str): string with information

    Raises:
        IndexError: if format is does not follow syntax.

    Returns:
        Metabolite: New Metabolite with given information.
    """
    line = [part.strip().rstrip() for part in line_string.split(",")]
    try:
        identifier = line[0]
        name = line[1]
        compartment = line[2]
        formula = line[3]
        charge = line[4]
    except IndexError:
        raise WrongSyntax(
            "Given line is invalid. Format is: id, name, compartment, "
            "chemical_formula, molecular_charge"
        )
    return _build_metabolite(
        identifier=identifier,
        name=name,
        compartment=compartment,
        charge=float(charge),
        formula=formula,
    )


# TODO: change it to identify if name should be formatted or reverse-formatted
def _fix_name(name: str) -> str:
    """
    Replaces hyphens in given name to underscores. Double hyphens are
    transformed into single underscores
    """
    return name.replace("-", "_")


def _get_metabolite(
    metabolite_dict: dict, compartment: str, model: Model = Model(0)
) -> Metabolite:
    """
    Return Metabolite based on data saved as a dictionary. The function
    will look in given model if compound is already in the model under
    another translation (name).

    Args:
        metabolite_dict (dict): dictionary with data of metabolite
        compartment (str): location of the metabolite.
        model (Model): Model to search for equivalents. Defaults to empty Model

    Returns:
        Metabolite: New object based on a dictionary.

    Raises:
        TypeError: If dictionary data does not represent a metabolite.
    """
    if not metabolite_dict["FORMULA"]:
        # TODO: change Error
        raise WrongDataError(
            "Given dictionary does not correspond to a Metabolite."
        )
    identifier = metabolite_dict["ENTRY"]
    # Try to obtain available translation
    with suppress(NoIntersectFound, KeyError):
        new_identifier = _first_item(
            first=model.metabolites,
            second=metabolite_dict["XREF"],
            revert=True,
        )
        # Only return from model if compartment is the same, otherwise
        # KeyError
        new_identifier = f"{_fix_name(name=new_identifier)}_{compartment}"
        metabolite = model.metabolites.get_by_id(new_identifier)
        debug_log.warning(
            f"Metabolite '{metabolite_dict['ENTRY']}' found in given model "
            f"under '{new_identifier}'"
        )
        return metabolite
    # Format if above fails
    identifier = f"{_fix_name(name=identifier)}_{compartment}"
    # Create object, add logging and then return it.
    return _build_metabolite(
        identifier=identifier,
        formula=metabolite_dict["FORMULA"],
        name=metabolite_dict["NAME"],
        charge=metabolite_dict["CHARGE"],
        compartment=compartment,
    )


def _convert_string_metabolite(line: str, model: Model, **kwargs):
    """
    Transform a string into a Metabolite and returns it. This object can be
    either custom or from a database. If the compound if found under a
    different name, this will be returned instead.

    Custom metabolite syntax:

    :code:`formatted_identifier, name, compartment, chemical_formula,
    molecular_charge`

    Metabolite from database:

    :code:`metabolite_identifier, compartment`

    Args:
        line (str): string with information of metabolite
        model (Model): model to test

     Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.

    Returns:
        Metabolite: New metabolite object
    """
    if line.count(",") > 1:
        # TODO: to gain perfomance, search for name and then create Metabolite
        new_metabolite = _metabolite_from_string(line_string=line)
    else:
        # Retrieve from databse
        segment = (part.strip().rstrip() for part in line.split(sep=","))
        # FIXME: include multiple databases=_identify_database()
        metabolite_dict = get_data(
            identifier=next(segment), debug_level=10, **kwargs
        )
        new_metabolite = _get_metabolite(
            metabolite_dict=metabolite_dict,
            compartment=next(segment),
            model=model,
        )
    return new_metabolite


def _get_file_metabolites(
    model: Model, filename: Path, **kwargs
) -> List[Metabolite]:
    """
    Return a list with Metabolites in file. If found in model under a different
    name, it will be included in the list instead.

    Custom metabolite syntax:

    :code:`formatted_identifier, name, compartment, chemical_formula,
    molecular_charge`

    Metabolite from database:

    :code:`metabolite_identifier, compartment`

    Args:
        model (Model): model to test
        filename (Path): location of the file with metabolites

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.

    Raises:
        TypeError: if model is invalid
        FileNotFoundError: if given file is not found
    """
    if not filename.exists():
        raise FileNotFoundError(
            f'Given file in "{str(filename)}" does not exists.'
        )
    # For each line, build and add metabolite. If a Metabolite is no properly
    # created, either raise an Error or use a default.
    with open(filename, "r") as f:
        new_metabolites = list()
        for line in _read_lines(f=f):
            new_metabolites.append(
                _convert_string_metabolite(line=line, model=model, **kwargs)
            )
    return new_metabolites


def _return_duplicate(data_dict: dict) -> bool:
    """
    Check for the duplicate in a dictionary with prefixes-
    """
    sequence = Counter([item[2:] for item in data_dict.keys()])
    return sequence.most_common(1)[0][0]


def _get_reaction(
    data_dict: dict,
    compartment: str,
    replacement: dict,
    model: Model,
    **kwargs,
) -> Reaction:
    """
    Creates a Reaction object from given dictionary with data. Location of the
    reactions can be set with the argument 'compartment'. The method will look
    in given model if the reaction and/or their corresponding metabolite are
    already in the model with other identifiers.

    Args:
        data_dict (dict): Dictionary with data of a Reaction.
        compartment (str): Locations of the reactions
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
        model (Model): Model to search for equivalents.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.

    Returns:
        Reaction: New reaction based on dictionary
    """
    identifier = data_dict["ENTRY"]
    # Try to obtain if information is available
    with suppress(NoIntersectFound, KeyError):
        new_identifier = _first_item(
            first=model.reactions, second=data_dict["XREF"], revert=True
        )
        new_identifier = f"{_fix_name(name=new_identifier)}_{compartment}"
        reaction = model.reactions.get_by_id(new_identifier)
        debug_log.warning(
            f"Reaction '{data_dict['ENTRY']}' found in given model "
            f"under '{new_identifier}'"
        )
        return reaction
    # Otherwise create from scratch
    reaction = Reaction(
        id=f'{_fix_name(name=data_dict["ENTRY"])}_{compartment}',
        name=data_dict["NAME"],
    )
    for identifier, coef in data_dict["EQUATION"].items():
        # Get rid of prefix r_ and l_
        identifier = identifier[2:]
        # First, replacement, since the identifier can be already in model
        with suppress(KeyError):
            identifier = replacement[identifier]
            debug_log.warning(
                f'Metabolite "{identifier}" replaced by '
                f'"{replacement[identifier]}".'
            )
        # TODO: check if this part is necesary
        # Retrieve data for metabolite
        try:
            # get metabolites from model if possible.
            metabolite = model.metabolites.get_by_id(identifier)
        except KeyError:
            metabolite = get_data(
                identifier=identifier, debug_level=10, **kwargs
            )
        # Checking if transport reaction
        if (
            data_dict["TRANSPORT"]
            and coef < 0
            and _return_duplicate(data_dict=data_dict["EQUATION"])
            == identifier
        ):
            # FIXME: temporary setting to extracellular
            metabolite = _get_metabolite(
                metabolite_dict=metabolite, compartment="e", model=model
            )
        else:
            # No transport
            metabolite = _get_metabolite(
                metabolite_dict=metabolite,
                compartment=compartment,
                model=model,
            )
        reaction.add_metabolites(metabolites_to_add={metabolite: coef})
        reaction.bounds = data_dict["BOUNDS"]
    return reaction


def _add_reactions_check(model: Model, reactions: List[Reaction]):
    """
    Check function that adds given Reactions to given model if it does not
    contain the reaction. It logs the skipped reactions.
    """
    for member in reactions:
        if member.id in [reaction.id for reaction in model.reactions]:
            debug_log.warning(
                f'Reaction "{member.id}" was found in given model. Skipping'
            )
            continue
        model.add_reactions([member])
        debug_log.info(f'Reaction "{member.id}" was added to model')


def _obtain_reaction(
    model: Model,
    identifier: str,
    directory: Path,
    database: str,
    compartment: str,
    replacement: dict,
):
    """
    Return Reaction object from local directory or given database. The method
    will look in given model if the reaction and/or their corresponding
    metabolites are already in the model under other names.

    Args:
        model (Model): Model to add reactions and search for equivalents.
        identifier (str): Original identifier of the reaction.
        directory (Path): Directory to search data.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.
        compartment (str): Location of the reaction.
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
    """
    # Obtain data
    data_dict = get_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=10,
    )
    # Transform it
    reaction = _get_reaction(
        data_dict=data_dict,
        compartment=compartment,
        directory=directory,
        database=database,
        replacement=replacement,
        model=model,
    )
    return reaction


def _dict_from_string(string_list: list) -> dict:
    """
    For given list of strings, creates a dictionary where keys are the
    identifiers of metabolites while values represent their corresponding
    coefficients

    Syntax:

    :code:`id_metabolite1: coefficient, id_metabolite2:coefficient ...
    id_metaboliteX: coefficient`

    Identifier has to end with an underscore and a compartment:

    E.g **`OXYGEN-MOLECULE_c: -1`**

    Args:
        string_list (list): List with strings with information about the
            metabolites

    Raises:
        WrongSyntax: If coefficient might be missing.

    Returns:
        dict: Dictionary with identifiers and coefficients
    """
    metabolites = dict()
    # Each member should have a pair as "WATER: -1"
    for member in string_list:
        pair = [x.strip().rstrip() for x in member.split(":")]
        name = pair[0]
        # Must have a len of two, otherwise something is missing
        try:
            value = float(pair[1])
            metabolites[name] = value
        except ValueError:
            raise WrongSyntax(f"Coefficient might be missing for {name}")
    return metabolites


def _reaction_from_string(
    line_string: str, directory: Path, database: str, model: Model = Model(0)
) -> Reaction:
    """
    Returns a custom reaction from given string, which includes the information
    of the Reaction and its metabolites. If metabolites are not in given model,
    they will be retrieved from the specified database. Function  will also
    search for translated-metabolites in the model.

    Syntax:

    :code:`reaction_identifier, reaction_name | metabolite_identifier1:
    coefficient,`
    :code:`metabolite_identifier2:coefficient, ..., metabolite_identifierX:
    coefficient`

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`OXYGEN-MOLECULE_c: -1`**

    Else, include a model to retrieve identifiers from it.

    Args:
        line_string (str): string with information
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.
        model (Model, Optional): A model to obtain metabolite objects from.
            Defaults to an empty Model.

    Raises:
        WrongSyntax: if identifier has a wrong format.

    Returns:
        Reaction: New reaction object.
    """
    segments = [x.strip().rstrip() for x in line_string.split("|")]
    rxn_id_name = [x.strip().rstrip() for x in segments[0].split(",")]
    try:  # no blanks
        string_list = [x.strip().rstrip() for x in segments[1].split(",")]
        meta_dict = _dict_from_string(string_list=string_list)
    except IndexError:  # wrong format
        raise WrongSyntax(
            f'No delimiter "|" found for {rxn_id_name[0].split(",")[0]}.'
        )
    # In case of regular of regular syntax
    if len(rxn_id_name) == 2:
        rxn_id, rxnName = rxn_id_name
    elif len(rxn_id_name) == 1 and rxn_id_name[0] != "":
        rxn_id = rxnName = rxn_id_name[0]
    else:  # blank space
        raise WrongSyntax(f"Wrong format for {segments}. ID not detected")
    # Create Base reaction and then fill it with its components.
    new_reaction = Reaction(id=rxn_id, name=rxnName)
    for identifier, coef in meta_dict.items():
        # Either get from model, or retrieve it.
        try:
            metabolite = model.metabolites.get_by_id(identifier)
        except KeyError:
            # _get_metabolite will also search for the metabolite under a
            # different name.
            compartment = identifier[-1:]
            identifier = identifier[:-2]
            # It is necessary to build the metabolite.
            data_dict = get_data(
                directory=directory,
                identifier=identifier,
                database=database,
                debug_level=10,
            )
            metabolite = _get_metabolite(
                metabolite_dict=data_dict, compartment=compartment, model=model
            )
        new_reaction.add_metabolites({metabolite: coef})
        debug_log.debug(
            f'Metabolite "{metabolite.id}" added to Reaction '
            f'"{new_reaction.id}".'
        )
    # FIXME: making all reversible
    new_reaction.bounds = (-1000, 1000)
    return new_reaction


def _convert_string_reaction(
    line: str,
    model: Model,
    directory: Path,
    database: str,
    replacement: dict = {},
) -> Reaction:
    """
    Returns a Reaction from string. It can be either custom, or the identifier
    for a reaction in a database. The function will search for the reactions
    and its corresponding metabolites under other names inside the model and
    return it, if necessary.

    Syntax:

    :code:`reaction_identifier, compartment`

    Form custom reactions:

    :code:`reaction_identifier, reaction_name | metabolite_identifier1:
    coefficient,`
    :code:`metabolite_identifier2:coefficient, ..., metabolite_identifierX:
    coefficient`

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`OXYGEN-MOLECULE_c: -1`**

    Args:
        line (str): String with custom reaction or identifier of reaction
        model (Model): Model to search for object if necessary.
        directory (Path): Path to directory, where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.
    """
    try:
        # Create custom reaction
        new_reaction = _reaction_from_string(
            line_string=line,
            directory=directory,
            database=database,
            model=model,
        )
    # If delimiter is not found, then it must be a reaction
    except WrongSyntax:
        # add reaction from root. Get only left part
        seqment = (part.strip().rstrip() for part in line.split(","))
        identifier = next(seqment)
        with suppress(KeyError):
            identifier = replacement[identifier]
        # TODO: identify database
        new_reaction = _obtain_reaction(
            model=model,
            identifier=identifier,
            directory=directory,
            database=database,
            compartment=next(seqment),
            replacement=replacement,
        )
    return new_reaction


def _get_file_reactions(
    model: Model, filename: Path, **kwargs
) -> List[Reaction]:
    """
    Returns list with reactions from file. All reactions can be either created
    manually or retrieved from a database. For each reactions, its always
    checks for mass balance.

    Custom reactions:

    :code:`reaction_identifier, reaction_name | metabolite_identifier1:
    coefficient,`
    :code:`metabolite_identifier2:coefficient, ..., metabolite_identifierX:
    coefficient`

    From database:

    :code:`original_identifier, compartment`

    Identifier has to end with an underscore and a compartment:

    E.g **`OXYGEN-MOLECULE_c: -1`**

    Args:
        model (Model): model to test
        filename (Path): location of the file with reaction information

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.

    Raises:
        FileNotFoundError: if file does not exists
    """
    # TODO: add mass balance check
    if not filename.exists():
        raise FileNotFoundError(f"Filename {filename.name} does not exist")
    with open(filename, "r") as f:
        new_reactions = list()
        for line in _read_lines(f=f):
            new_reactions.append(
                _convert_string_reaction(line=line, model=model, **kwargs)
            )
    return new_reactions


def _ident_pathway(data_dict: dict) -> Generator:
    """
    Checks if given dictionary belongs to a pathway object. Returns its as
    generator if True. Else, raises WrongDataError.
    """
    if data_dict["TYPE"] == "Pathway":
        debug_log.info(
            f"Object '{data_dict['ENTRY']}' identified as a pathway"
        )
        yield data_dict
    else:
        raise WrongDataError("Data does not belong to a pathway")


def _ident_reaction(
    data_dict: dict,
    directory: Path,
    replacement: dict,
    database: str,
    compartment: str,
    show_imbalance: bool,
    stop_imbalance: bool,
    model: Model,
) -> Generator:
    """
    Tries to identify given dictionary if it includes information of a
    reaction. If True, it will yield a Reaction.

    Args:
        data_dict (dict): dictionary to examine.
        directory (Path): directory to retrived and store data.
        replacement (dict): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.
        database (str): name of the database to query retrive data.
        compartment: location of the reaction.
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
            Defaults to False.
        show_imbalance (bool): If unbalanced reactionis found, show output.
            Defaults to True.

    Returns:
       Generator: Reaction from dictionary.

    Raises:
        WrongDataError: If data does not include reaction information.
    """
    try:
        # First create object and then yield it. Otherwise, object will not be
        # created correctly
        reaction = _get_reaction(
            data_dict=data_dict,
            directory=directory,
            replacement=replacement,
            database=database,
            compartment=compartment,
            model=model,
        )
        check_imbalance(
            reaction=reaction,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
        )
        debug_log.info(
            f"Object '{data_dict['ENTRY']}' identified as a reaction"
        )
        yield reaction
    except KeyError:
        raise WrongDataError("Data does not belong to a reaction.")


def _ident_metabolite(
    data_dict: dict, compartment: str, model: Model
) -> Generator:
    """
    Tries to identify given dictionary if its includes metabolite information,
    and transforms it into a metabolite in a generator

    Args:
        data_dict: dictionary with data of a metabolite.
        compartment: location of the metabolite.

    Returns:
        Generator: new Metabolite based on given dictionary.

    Raises:
        WrongDataError: If data does not include information of a metabolite.
    """
    try:
        # First create object and then yield it. Otherwise, object will not be
        # created correctly
        metabolite = _get_metabolite(
            metabolite_dict=data_dict, compartment=compartment, model=model
        )
        debug_log.info(f"Object '{metabolite.id}' identified as a metabolite")
        yield metabolite
    except KeyError:
        raise WrongDataError("Data does not belong to a metabolite")


def create_object(
    identifier: str,
    directory: Path,
    database: str,
    compartment: str,
    replacement: dict = {},
    show_imbalance: bool = True,
    stop_imbalance: bool = False,
    model: Model = Model(0),
    model_id="universal",
) -> Union[Reaction, Metabolite, dict]:
    """
    Creates and returns COBRApy object based on given identifier and database.
    Identifier names will be formatted.

    .. hint:: Hyphens will become underscores. Double hyphens become single\
    underscores.

    Args:
        identifier (str): Original identifier for database
        directory (Path): Path to directory where data is stored.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.
        compartment (str): Location of the object. In case of reaction, all
            metabolites will be included in the same location.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}. Does not apply to
            pathways.

    Arguments for reactions:
        stop_imbalance (bool, optional): If unbalanced reaction is found, stop
            process. Defaults to False.
        show_imbalance (bool, optional): If unbalanced reaction is found, show
            output. Defaults to True.
        model (Model, optional): Model to add search for translated metabolites
            or reactions. Defaults to a empty model.
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathway are not available.
            Defaults to: "universal"

    Returns:
        Union[Reaction, Metabolite]: A Reaction or Metabolite object; or the
            information for a pathway.
    """
    data_dict = get_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=10,
        model_id=model_id,
    )
    # Since it is only a single item, next() can be used
    for method in (
        _ident_pathway(data_dict=data_dict),
        _ident_metabolite(
            data_dict=data_dict, compartment=compartment, model=model
        ),
        _ident_reaction(
            data_dict=data_dict,
            directory=directory,
            replacement=replacement,
            database=database,
            compartment=compartment,
            show_imbalance=show_imbalance,
            stop_imbalance=stop_imbalance,
            model=model,
        ),
    ):
        # Try to return object, unless it cannot be identified
        with suppress(WrongDataError):
            return next(method)
    raise Warning("Data cannot be identified. Examine with 'get_data'.")


def __add_metabolites_check(model: Model, metabolites: List[Metabolite]):
    """
    Checks if given metabolites are already in the model. If not, they will be
    added into the model.

    Args:
        model (Model): Model to extend.
        metabolites (List[Metabolites]): List with Metabolites.
    """
    # A Loop in necessary to log the skipped metabolites.
    for member in metabolites:
        if member.id not in [
            metabolite.id for metabolite in model.metabolites
        ]:
            model.add_metabolites(metabolite_list=member)
            debug_log.info(f'Metabolite "{member.id}" was added to model.')
            continue
        debug_log.warning(
            f'Metabolite "{member.id}" was already in model. Skipping.'
        )


def add_metabolites(model: Model, obj: Any, **kwargs):
    """
    Adds given object into the model. The options are:

     - Path: A file with components. E. g:
        Path.cwd().joinpath("file_with_names.txt")
     - Metabolite: A single Metabolite.
     - List[Metabolites]: A list with multiple Metabolite objects.
     - str: Either the identifier with its corresponding compartment or a
     string with the whole attributes. This applies for the Path option. E.g:

        Custom metabolite syntax:

        :code:`formatted_identifier, name, compartment, chemical_formula,
        molecular_charge`

        Metabolite from database:

        :code:`metabolite_identifier, compartment`

     - List[str]: A list with multiple str with the mentioned syntax.

    Args:
        model (Model): Model to expand and search for metabolites.
        obj: A Path; a list with either strings or Metabolite objects,
            or a single string. See syntax above.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.

    Raises:
        WrongSyntax (from str): If syntax is not followed correctly as
            mentioned above.
        ValueError: If Keyword Arguments are missing.
        FileNotFoundError (from Path): if file does not exists
    """
    try:
        # In case of a Path
        if isinstance(obj, Path):
            # These variable will raise KeyError if no kwargs are passed.
            directory = kwargs["directory"]
            database = kwargs["database"]
            new_metabolites = _get_file_metabolites(
                # as Path
                filename=obj,
                model=model,
                directory=directory,
                database=database,
            )
        # In case of a single Metabolite
        elif isinstance(obj, Metabolite):
            new_metabolites = [obj]
        # Unless, iterable with Metabolites
        elif all([isinstance(member, Metabolite) for member in obj]):
            new_metabolites = obj
        # or a list with str
        elif all([isinstance(member, str) for member in obj]) or isinstance(
            obj, str
        ):
            # These variable will raise KeyError if no kwargs are passed.
            directory = kwargs["directory"]
            database = kwargs["database"]
            # Make a list
            if isinstance(obj, str):
                obj = [obj]
            # Create new list
            new_metabolites = [
                _convert_string_metabolite(
                    model=model,
                    line=line,
                    directory=directory,
                    database=database,
                )
                for line in obj
            ]
        # Raise error if wrong
        else:
            raise WrongDataError(
                "Given object is not the Type mentioned in the docstrings."
            )
        # Otherwise, it must be a list with Metabolites.
        __add_metabolites_check(model=model, metabolites=new_metabolites)
    except KeyError:
        raise ValueError("Keyword Arguments are missing for given object.")


def __add_reactions_check(model: Model, reactions: List[Reaction]):
    """
    Checks if given reactions are already in the model. If not, they will be
    added into the model.

    Args:
        model (Model): Model to extend.
        reactions (List[Metabolites]): List with Reactions.
    """
    # A Loop in necessary to log the skipped metabolites.
    for member in reactions:
        if member.id not in [reaction.id for reaction in model.reactions]:
            model.add_reactions(reaction_list=[member])
            debug_log.info(f'Reaction "{member.id}" was added to model.')
            continue
        debug_log.warning(
            f'Reaction "{member.id}" was already in model. Skipping.'
        )


def add_reactions(model: Model, obj: Any, **kwargs):
    """Adds given object into the model. The options are:

     - Path: A file with components. E. g:
        Path.cwd().joinpath("file_with_names.txt")
     - List[Reactions]: A list with regular Reactions
     - str: Either the identifier with its corresponding compartment or a
     string with the whole components. This appplies for the Path option. E.g:

        :code:`reaction_identifier, compartment`

        For custom reactions

        :code:`reaction_identifier, reaction_name | metabolite_identifier1:
        coefficient,`
        :code:`metabolite_identifier2:coefficient, ..., metabolite_identifierX:
        coefficient`

        Identifiers of metabolites have to end with an underscore and a
        compartment:

            E.g  **`OXYGEN-MOLECULE_c: -1`**

     - List[str]: A list with multiple str with the mentioned syntax.

    Args:
        model (Model): Model to expand and search for reactions.
        obj: A Path; a list with either strings or Reaction objects,
            or a single string. See syntax above.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.
        replacement (dict): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.

    Raises:
        WrongSyntax (from str): If syntax is not followed correctly as
            mentioned above.
        ValueError: If Keyword Arguments are missing.
        FileNotFoundError (from Path): if file does not exists
    """
    try:
        # In case of a Path
        if isinstance(obj, Path):
            # These variable will raise KeyError if no kwargs are passed.
            directory = kwargs["directory"]
            database = kwargs["database"]
            # Empty dictionary is nothing assigned
            replacement = kwargs.pop("replacement", dict())
            new_reactions = _get_file_reactions(
                model=model,
                filename=obj,
                directory=directory,
                database=database,
                replacement=replacement,
            )
        # In case of single Reaction
        elif isinstance(obj, Reaction):
            new_reactions = [obj]
        # Unless, iterable with Reactions.
        elif all([isinstance(member, Reaction) for member in obj]):
            new_reactions = obj
        # or a list with str
        elif all([isinstance(member, str) for member in obj]) or isinstance(
            obj, str
        ):
            directory = kwargs["directory"]
            database = kwargs["database"]
            # Empty dictionary is nothing assigned
            replacement = kwargs.pop("replacement", dict())
            # Make a list
            if isinstance(obj, str):
                obj = [obj]
            new_reactions = [
                _convert_string_reaction(
                    line=line,
                    model=model,
                    directory=directory,
                    database=database,
                    replacement=replacement,
                )
                for line in obj
            ]
        # Raise error if wrong
        else:
            raise WrongDataError(
                "Given object is not the Type mentioned in the docstrings."
            )
        # Otherwise, it must be a list with Metabolites.
        __add_reactions_check(model=model, reactions=new_reactions)
    except KeyError:
        raise ValueError("Keyword Arguments are missing for given object.")
