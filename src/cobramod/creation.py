#!/usr/bin/env python3
"""Object creation

This module handles the creation of COBRApy's objects
:func:`cobra.core.metabolite.Metabolite` and
:func:`cobra.core.reaction.Reaction`. Dictionaries with the data from given
database is used. Important functions are:

- create_object: Creates and Returns a COBRApy object.
- add_reactions_from_file: Add reactions in a file to a model.
- _custom_reaction: Return a custom reaction
- add_reaction:
- _include_file_metabolite: Add metabolites in a file to a model.
"""
from collections import Counter
from contextlib import suppress
from pathlib import Path
from typing import Union, Generator

from cobra import Metabolite, Model, Reaction

from cobramod.debug import debug_log
from cobramod.error import WrongDataError, NoIntersectFound
from cobramod.mod_parser import get_data
from cobramod.utils import _read_lines, check_imbalance, _first_item


def _metabolite_from_string(line_string: str) -> Metabolite:
    """
    Creates a Metabolite object based on a string.

    The string must follow the syntax:

    :code:`formatted identifier, name , compartment, chemical_formula,
    molecular_charge`

    Args:
        line_string (str): string with information

    Raises:
        TypeError: if no string is identifier
        IndexError: if format is invalid

    Returns:
        Metabolite: New Metabolite with given information.
    """
    if not isinstance(line_string, str):
        raise TypeError("Argument must be a str")
    line = [part.strip().rstrip() for part in line_string.split(",")]
    try:
        meta_id = line[0]
        meta_name = line[1]
        meta_comp = line[2]
        meta_formula = line[3]
        meta_charge = line[4]
    except IndexError:
        raise IndexError(
            "Given line is invalid. Format is: id, name, compartment, "
            "chemical_formula, molecular_charge"
        )
    return Metabolite(
        id=meta_id,
        name=meta_name,
        compartment=meta_comp,
        charge=meta_charge,
        formula=meta_formula,
    )


# TODO: change it to identify if name should be formatted or reverse-formatted
def _fix_name(name: str) -> str:
    """
    Replaces hyphens in given name to underscores. Double hyphens are
    transformed into single underscores
    """
    return name.replace("-", "_")


def _build_metabolite(
    metabolite_dict: dict, compartment: str, model: Model = Model(0)
) -> Metabolite:
    """
    Builds and return a metabolite based on data saved as a dictionary. The
    method will look in given model if compound is already in the model with
    another translation.

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
        raise TypeError(
            "Given dictionary does not correspond to a Metabolite."
        )
    identifier = metabolite_dict["ENTRY"]
    # Try to obtain available translation
    with suppress(NoIntersectFound, KeyError):
        identifier = _first_item(
            first=model.metabolites,
            second=metabolite_dict["XREF"],
            revert=True,
        )
        # Only return from model if compartment is the same, otherwise
        # KeyError
        identifier = f"{_fix_name(name=identifier)}_{compartment}"
        debug_log.warning(
            f"Metabolite '{metabolite_dict['ENTRY']}' found in given model "
            f"under '{identifier}'"
        )
        metabolite = model.metabolites.get_by_id(identifier)
        return metabolite
    # Format if above fails
    identifier = f"{_fix_name(name=identifier)}_{compartment}"
    debug_log.debug(f'Metabolite "{identifier}" created.')
    return Metabolite(
        id=identifier,
        formula=metabolite_dict["FORMULA"],
        name=metabolite_dict["NAME"],
        charge=metabolite_dict["CHARGE"],
        compartment=compartment,
    )


def _check_if_meta_in_model(metabolite: str, model: Model) -> bool:
    """
    Returns if metabolite identifier is found in given model.
    """
    return metabolite in [meta.id for meta in model.metabolites]


def _add_if_not_found_model(model: Model, metabolite: Metabolite):
    """
    Checks if given Metabolite object is found in given Model. If not, it will
    be added

    Args:
        model (Model): model to test
        metabolite (Metabolite): Metabolite object to test
    """
    if _check_if_meta_in_model(model=model, metabolite=metabolite.id):
        debug_log.warning(
            f'Metabolite "{metabolite.id}" was found in given model. Skipping'
        )
    else:
        model.add_metabolites([metabolite])
        debug_log.info(f'Metabolite "{metabolite.id}" was added to model')


def _include_string(line: str, model: Model, **kwargs):
    """
    Transform a string into a Metabolite object and appends it into model.
    The Metabolite can be either custom or from a database. Returns new
    Metabolite object.

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
        seqment = (part.strip().rstrip() for part in line.split(sep=","))
        # FIXME: include multiple databases=_identify_database()
        metabolite_dict = get_data(
            identifier=next(seqment), debug_level=10, **kwargs
        )
        new_metabolite = _build_metabolite(
            metabolite_dict=metabolite_dict,
            compartment=next(seqment),
            model=model,
        )
    _add_if_not_found_model(model=model, metabolite=new_metabolite)
    return new_metabolite


def _include_file_metabolite(model: Model, filename: Path, **kwargs):
    """
    Creates and add new Metabolites specified in given file to given model.

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
        raise FileNotFoundError
    # For each line, build and add metabolite. If a Metabolite is no properly
    # created, either raise an Error or use a default.
    with open(filename, "r") as f:
        for line in _read_lines(f=f):
            try:
                _include_string(line=line, model=model, **kwargs)
            # TODO: create new Error
            except Warning:
                # FIXME: add function to search for common missing databases
                kwargs["database"] = "META"
                _include_string(line=line, model=model, **kwargs)


def _return_duplicate(data_dict: dict) -> bool:
    """
    Check for the duplicate in a dictionary with prefixes-
    """
    sequence = Counter([item[2:] for item in data_dict.keys()])
    return sequence.most_common(1)[0][0]


def _build_reaction(
    data_dict: dict,
    compartment: str,
    replacement: dict,
    model: Model,
    **kwargs,
) -> Reaction:
    """
    Creates a Reaction object from given dictionary with data. Location of the
    reactions can be set with the argument 'compartment'.  Metabolites can be
    replaced by using the dictionary 'replacement' with the following
    syntax:

    :code:`old_identifier: replacement`

    The method will look in given model if the reaction and/or their
    corresponding metabolite are already in the model with other identifiers.

    Args:
        data_dict (dict): Dictionary with data of a Reaction.
        compartment (str): Locations of the reactions
        replacement (dict): original identifiers to be replaced.
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
        identifier = _first_item(
            first=model.reactions, second=data_dict["XREF"], revert=True
        )
        identifier = f"{_fix_name(name=identifier)}_{compartment}"
        debug_log.warning(
            f"Reaction '{data_dict['ENTRY']}' found in given model "
            f"under '{identifier}'"
        )
        return model.reactions.get_by_id(identifier)
    # Otherwise create from scratch
    reaction = Reaction(
        id=f'{_fix_name(name=data_dict["ENTRY"])}_{compartment}',
        name=data_dict["NAME"],
    )
    for identifier, coef in data_dict["EQUATION"].items():
        # Get rid of prefix r_ and l_
        identifier = identifier[2:]
        # TODO: add option to get metabolites from model
        # First, replacement, since the identifier can be already in model
        with suppress(KeyError):
            identifier = replacement[identifier]
            debug_log.warning(
                f'Metabolite "{identifier}" replaced by '
                f'"{replacement[identifier]}".'
            )
        # TODO: check if this part is necesary
        # Retrieve data for metabolite
        metabolite = get_data(identifier=identifier, debug_level=10, **kwargs)
        # Checking if transport reaction
        if (
            data_dict["TRANSPORT"]
            and coef < 0
            and _return_duplicate(data_dict=data_dict["EQUATION"])
            == identifier
        ):
            # FIXME: temporary setting to extracellular
            metabolite = _build_metabolite(
                metabolite_dict=metabolite, compartment="e", model=model
            )
        else:
            # No transport
            metabolite = _build_metabolite(
                metabolite_dict=metabolite,
                compartment=compartment,
                model=model,
            )
        reaction.add_metabolites(metabolites_to_add={metabolite: coef})
        reaction.bounds = data_dict["BOUNDS"]
    return reaction


def _add_reaction_check(model: Model, reaction: Reaction):
    """
    Check function that adds given Reaction to given model if it does not
    contain the reaction. Function checks for identifier of reaction.
    """
    if reaction.id in [reaction.id for reaction in model.reactions]:
        debug_log.warning(
            f'Reaction "{reaction.id}" was found in given model. Skipping'
        )
        return
    model.add_reactions([reaction])
    debug_log.info(f'Reaction "{reaction.id}" was added to model')


def _include_reaction(
    model: Model,
    identifier: str,
    directory: Path,
    database: str,
    compartment: str,
    replacement: dict,
):
    """
    Creates and adds a Reaction object based on given identifier from given
    database. Metabolites can be replaced by using the dictionary
    replacement with the following syntax:

    :code:`old_identifier: replacement`

    The method will look in given model if the reaction and/or their
    corresponding metabolite are already in the model with other identifiers.

    Args:
        model (Model): Model to add reactions and search for equivalents.
        identifier (str): Original identifier of the reaction.
        directory (Path): Directory to search data.
        database (str): Name of database. Check
            :func:`cobramod.available_databases` for a list of names.
        compartment (str): Location of the reaction.
        replacement (dict): Original identifiers to be replaced.
    """
    # Obtain data
    data_dict = get_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=10,
    )
    # Transform it
    reaction = _build_reaction(
        data_dict=data_dict,
        compartment=compartment,
        directory=directory,
        database=database,
        replacement=replacement,
        model=model,
    )
    # Add to model
    # TODO: check if this can be restructured
    _add_reaction_check(model=model, reaction=reaction)


def _has_delimiter(line_string: str) -> bool:
    """
    Returns true if given string has a vertical bar '|'
    """
    return "|" in line_string


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
        TypeError: if format is wrong
        ValueError: if coefficient is missing

    Returns:
        dict: Dictionary with identifiers and coefficients
    """
    if not isinstance(string_list, list):
        raise TypeError("Line format is wrong")
    meta_dict = dict()
    for single in string_list:
        single = [x.strip().rstrip() for x in single.split(":")]
        meta_name = single[0]
        try:
            meta_value = float(single[1])
            meta_dict[meta_name] = meta_value
        except ValueError:
            raise ValueError(f"Coefficient might be missing for {meta_name}")
    return meta_dict


def _custom_reaction(
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
        IndexError: if not  '|' is not found.
        Warning: if identifier has a wrong format.

    Returns:
        Reaction: New reaction object.
    """
    segments = [x.strip().rstrip() for x in line_string.split("|")]
    rxn_id_name = [x.strip().rstrip() for x in segments[0].split(",")]
    try:  # no blanks
        string_list = [x.strip().rstrip() for x in segments[1].split(",")]
        meta_dict = _dict_from_string(string_list=string_list)
    except IndexError:  # wrong format
        # TODO: Create proper error
        raise IndexError(
            f'No delimiter "|" found for {rxn_id_name[0].split(",")[0]}.'
        )
    # In case of regular of regular syntax
    if len(rxn_id_name) == 2:
        rxn_id, rxnName = rxn_id_name
    # Otherwise, share same name and identifier
    elif len(rxn_id_name) == 1 and rxn_id_name[0] != "":
        rxn_id = rxnName = rxn_id_name[0]
    else:  # blank space
        # TODO: Create proper error
        raise Warning(f"Wrong format for {segments}. ID not detected")
    # Create Base reaction and then fill it with its components.
    new_reaction = Reaction(id=rxn_id, name=rxnName)
    for identifier, coef in meta_dict.items():
        try:
            # To avoid creation of metabolite if already in model
            metabolite = model.metabolites.get_by_id(identifier)
        except KeyError:
            compartment = identifier[-1:]
            identifier = identifier[:-2]
            data_dict = get_data(
                directory=directory,
                identifier=identifier,
                database=database,
                debug_level=10,
            )
            metabolite = _build_metabolite(
                metabolite_dict=data_dict, compartment=compartment, model=model
            )
        new_reaction.add_metabolites({metabolite: coef})
    # FIXME: making all reversible
    new_reaction.bounds = (-1000, 1000)
    return new_reaction


# TODO: This will be deprecated. It will become add_reactions
def _add_reaction_line_to_model(
    line: str,
    model: Model,
    directory: Path,
    database: str,
    replacement: dict = {},
    **kwargs,
):
    """
    From given string, it will identify if a custom Reaction or a Reaction
    from root can be created. It will build the reaction and adds it to
    given model

    Args:
        line (str): string with information
        model (Model): model to test
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA". "BIGG",
            "KEGG"

    Keyword Arguments:
        **kwargs: same as in :func:`cobramod.creation._custom_reaction`
    """
    if _has_delimiter(line_string=line):
        # create custom reaction
        new_reaction = _custom_reaction(
            line_string=line,
            directory=directory,
            database=database,
            model=model,
        )
        _add_reaction_check(model=model, reaction=new_reaction)
    else:
        # add reaction from root. Get only left part
        seqment = (part.strip().rstrip() for part in line.split(","))
        identifier = next(seqment)
        with suppress(KeyError):
            identifier = replacement[identifier]
        # TODO: identify database
        _include_reaction(
            model=model,
            identifier=identifier,
            directory=directory,
            database=database,
            compartment=next(seqment),
            replacement=replacement,
        )


def add_reactions_from_file(model: Model, filename: Path, **kwargs):
    """
    Adds new reactions to given Model. All reactions can be either created
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
        database (str): Name of database. Options: "META", "ARA".
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.

    Raises:
        FileNotFoundError: is file does not exists
    """
    # TODO: add mass balance check
    if not filename.exists():
        raise FileNotFoundError(f"Filename {filename.name} does not exist")
    with open(filename, "r") as f:
        for line in _read_lines(f=f):
            try:
                _add_reaction_line_to_model(line=line, model=model, **kwargs)
            except Warning:
                # TODO: add test
                kwargs["database"] = "META"
                _add_reaction_line_to_model(line=line, model=model, **kwargs)


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
        reaction = _build_reaction(
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
        metabolite = _build_metabolite(
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
        database (str): Name of the database. Options are: "META", "ARA",
            "KEGG", "BIGG"
        compartment (str): Location of the object. In case of reaction, all
            metabolites will be included in the same location.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}. Does not apply to
            pathways.
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


def add_reactions():
    pass
