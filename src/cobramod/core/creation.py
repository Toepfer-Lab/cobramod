#!/usr/bin/env python3
"""Object creation

This module handles the creation of COBRApy's objects
:class:`cobra.core.metabolite.Metabolite` and
:class:`cobra.core.reaction.Reaction`. Dictionaries containing the data of a
given database are used. Important functions are:

- create_object: Creates and Returns a COBRApy object.
- add_reactions: Add reactions from multiple sources.
- add_metabolites: Add reactions from multiple sources.

These functions are a mix of multiple simpler functions:

- _metabolite_from_string, _reaction_from_string: create objects from strings.
- _get_metabolite, _get_reaction: create objects from dictionary.
- _convert_string_reaction, _convert_string_metabolite: create objects from
  files.
"""
from collections import Counter
from contextlib import suppress
from pathlib import Path
from typing import Union, Generator, List, Any, Optional
from warnings import warn

from cobra import Metabolite, Model, Reaction
from requests import HTTPError

from cobramod.debug import debug_log
from cobramod.error import (
    WrongDataError,
    NoIntersectFound,
    WrongSyntax,
    WrongParserError,
    NoDelimiter,
)
from cobramod.core.retrieval import get_data
from cobramod.utils import _read_lines, check_imbalance, _first_item
from cobramod.core.genes import _genes_to_reaction


def _build_metabolite(
    identifier: str, formula: str, name: str, charge: float, compartment: str
):
    """
    Returns a basic :class:`cobra.Metabolite`. It will log a with a DEBUG
    level.

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
    debug_log.debug(f'Manually curated metabolite "{identifier}" was created.')
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
        debug_log.debug(
            f'Metabolite "{identifier}" was identified as a manually. '
            + "curated metabolite."
        )
    except IndexError:
        raise WrongSyntax(
            "Given line is invalid. Format is:\nid, name, compartment, "
            f"chemical_formula, molecular_charge.\nVerify line:\n{line_string}"
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
    Returns new string, where hyphens are replaced to underscores
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
        WrongDataError: If dictionary data does not represent a metabolite.
    """
    if not metabolite_dict["FORMULA"]:
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
        msg = (
            f"Metabolite '{metabolite_dict['ENTRY']}' was found as "
            f"'{new_identifier}'. Please curate if necessary."
        )
        debug_log.warning(msg=msg)
        warn(message=msg, category=UserWarning)
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


def _convert_string_metabolite(
    line: str, model: Model, **kwargs
) -> Metabolite:
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
            :obj:`cobramod.available_databases` for a list of names.

    Returns:
        Metabolite: New metabolite object
    """
    if line.count(",") > 1:
        # TODO: to gain performance, search for name and then create Metabolite
        new_metabolite = _metabolite_from_string(line_string=line)
    else:
        # Retrieve from database
        segment = (part.strip().rstrip() for part in line.split(sep=","))
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
    Return a list with Metabolites from a file. If a metabolite is found in
    model under a different name, it will be included in the list instead.

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
            :obj:`cobramod.available_databases` for a list of names.

    Raises:
        TypeError: if model is invalid
        FileNotFoundError: if given file is not found
    """
    if not filename.exists():
        raise FileNotFoundError(
            f'Given file in "{str(filename)}" does not exists. '
            + "Please create the given file."
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
    show_imbalance: bool,
    stop_imbalance: bool,
    **kwargs,
) -> Reaction:
    """
    Creates a Reaction object from given dictionary with data. Location of the
    reactions can be set with the argument 'compartment'. The method will look
    in given model if the reaction and/or their corresponding metabolite are
    already in the model with other identifiers. This method will create the
    corresponding genes and add the to the reaction.

    Args:
        data_dict (dict): Dictionary with data of a Reaction.
        compartment (str): Locations of the reactions
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
        model (Model): Model to search for equivalents.
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names.

    Returns:
        Reaction: New reaction based on dictionary
    """
    identifier = data_dict["ENTRY"]
    if data_dict["DATABASE"] == "META":
        msg = (
            f'Metabolic pathway information for reaction "{identifier}" comes '
            'from database "META". No genes will be created. Please use a '
            "specie-specific sub-database from BioCyc."
        )
        debug_log.warning(msg)
        warn(message=msg, category=UserWarning)

    # Try to obtain if information is available
    with suppress(NoIntersectFound, KeyError):
        new_identifier = _first_item(
            first=model.reactions, second=data_dict["XREF"], revert=True
        )
        new_identifier = f"{_fix_name(name=new_identifier)}_{compartment}"
        reaction = model.reactions.get_by_id(new_identifier)
        msg = (
            f"Reaction '{data_dict['ENTRY']}' was found as "
            f"'{new_identifier}'. Please curate if necessary"
        )
        debug_log.warning(msg=msg)
        warn(message=msg, category=UserWarning)
        return reaction
    # Otherwise create from scratch
    reaction = Reaction(
        id=f'{_fix_name(name=data_dict["ENTRY"])}_{compartment}',
        name=data_dict["NAME"],
    )
    for identifier, coefficient in data_dict["EQUATION"].items():
        # Get rid of prefix r_ and l_
        identifier = identifier[2:]
        # First, replacement, since the identifier can be already in model
        with suppress(KeyError):
            identifier = replacement[identifier]
            msg = (
                f'Metabolite "{identifier}" found in "replacement" '
                f"Metabolite will be replaced by "
                f'"{replacement[identifier]}" in reaction "{reaction.id}".'
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
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
            and coefficient < 0
            and _return_duplicate(data_dict=data_dict["EQUATION"])
            == identifier
        ):
            # Temporary setting to extracellular
            metabolite = _get_metabolite(
                metabolite_dict=metabolite, compartment="e", model=model
            )
            msg = (
                f'Reaction "{reaction.id}" has metabolite "{metabolite.id}" on'
                " both sides of the equation (e.g transport reaction). COBRApy"
                " ignores these metabolites. To avoid this, by default, "
                "CobraMod will assign one metabolite to the extracellular "
                "compartment. Please curate the reaction if necessary."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
        else:
            # No transport
            metabolite = _get_metabolite(
                metabolite_dict=metabolite,
                compartment=compartment,
                model=model,
            )
        reaction.add_metabolites(metabolites_to_add={metabolite: coefficient})
        reaction.bounds = data_dict["BOUNDS"]
    # Add genes and check imbalance
    _genes_to_reaction(reaction=reaction, data_dict=data_dict)
    check_imbalance(
        reaction=reaction,
        show_imbalance=show_imbalance,
        stop_imbalance=stop_imbalance,
    )
    return reaction


def _add_reactions_check(model: Model, reactions: List[Reaction]):
    """
    Check function that adds given Reactions to given model if it does not
    contain the reaction. It logs the skipped reactions.
    """
    for member in reactions:
        if member.id in [reaction.id for reaction in model.reactions]:
            msg = (
                f'Reaction "{member.id}" is already present in the model. '
                "Skipping additions."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
            continue
        model.add_reactions([member])
        debug_log.info(f'Reaction "{member.id}" was added to model.')


def _obtain_reaction(
    model: Model,
    identifier: str,
    directory: Path,
    database: str,
    compartment: str,
    replacement: dict,
    stop_imbalance: bool,
    show_imbalance: bool,
    genome: Optional[str],
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
            :obj:`cobramod.available_databases` for a list of names.
        compartment (str): Location of the reaction.
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers.
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        genome (str): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie.
            List available at https://www.genome.jp/kegg/catalog/org_list.html
    """
    # Obtain data
    data_dict = get_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=10,
        genome=genome,
    )
    # Transform it
    reaction = _get_reaction(
        data_dict=data_dict,
        compartment=compartment,
        directory=directory,
        database=database,
        replacement=replacement,
        model=model,
        show_imbalance=show_imbalance,
        stop_imbalance=stop_imbalance,
    )
    return reaction


def _reaction_information(string: str) -> tuple:
    """
    Parses the string and transforms it to dictionary with the metabolite and
    their corresponding coefficients. It will raise WrongSyntax if the format
    is wrong. Returns the tuple with the bounds of the reactions and its
    participants as a dictionary. Reversibility will depend on the type of
    arrow that. Options are:

    **`<--, -->, <=>, <->`**
    """
    # Define bounds
    arrows = {
        "-->": (0, 1000),
        "<--": (-1000, 0),
        "<=>": (-1000, 1000),
        "<->": (-1000, 1000),
    }
    for separator in arrows.keys():
        if separator in string:
            bounds = arrows[separator]
            break
    # Verify bounds
    try:
        bounds
    except UnboundLocalError:
        raise WrongSyntax(
            f"Given string does not have a proper arrow.\n{string}\n"
            f'Please use one of this options:\n"{arrows.keys()}"'
        )
    # Separate parts
    reactants, products = [
        part.rstrip().strip() for part in string.split(separator)
    ]
    participants = dict()
    # Check for each side and make a dictionary with the coefficients
    for side in (reactants, products):
        FACTOR = 1
        if side is reactants:
            FACTOR = -1
        for metabolite in side.split("+"):
            # Get rid off blanks in the sides
            metabolite = metabolite.rstrip().strip()
            # Obtain coefficient or give a 1 to the metabolite if nothing found
            try:
                coefficient, identifier = [
                    string.strip().rstrip() for string in metabolite.split(" ")
                ]
                if "_" not in identifier[-2:]:
                    raise WrongSyntax(
                        f'Metabolite "{identifier}" is missing a compartment. '
                        + "Please add an underscore and the corresponding"
                        + "compartment. E.g. WATER_p"
                    )
            except ValueError:
                coefficient = "1"
                identifier = metabolite
            try:
                # Make sure that no empty strings are added.
                if identifier:
                    participants[identifier] = float(coefficient) * FACTOR
            except ValueError:
                raise WrongSyntax(
                    f"Check the syntax for the given line:\n{string}. "
                    + "Please use COBRApy standards. Example:\n"
                    + "g6p_c --> f6p_c + h2o_c"
                )
    return bounds, participants


def _reaction_from_string(
    line_string: str,
    directory: Path,
    database: str,
    stop_imbalance: bool,
    show_imbalance: bool,
    model: Model = Model(0),
) -> Reaction:
    """
    Returns a custom reaction from given string, which includes the information
    of the Reaction and its metabolites. If metabolites are not in given model,
    they will be retrieved from the specified database. Function  will also
    search for translated-metabolites in the model.

    Syntax:

    :code:`reaction_identifier, reaction_name | coefficient metabolite <->
    coefficient metabolite

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`4 OXYGEN-MOLECULE_c`**

    Else, include a model to retrieve identifiers from it.

    Reversibility will depend on the type of arrow that. Options are

    **`<--, -->, <=>, <->`**

    Args:
        line_string (str): string with information
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names.
        model (Model, Optional): A model to obtain metabolite objects from.
            Defaults to an empty Model.
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.

    Raises:
        WrongSyntax: if identifier has a wrong format.
        NoDelimiter: If no delimiter "|" is found.

    Returns:
        Reaction: New reaction object.
    """
    if "|" not in line_string:
        raise NoDelimiter(string=line_string)
    # Remove blank and obtains information
    info, reaction = [
        string.strip().rstrip() for string in line_string.split("|")
    ]
    # It is possible to pass only the identifier
    if "," in info:
        reaction_id, reaction_name = [
            string.strip().rstrip() for string in info.split(",")
        ]
    else:
        reaction_id = reaction_name = info
    # Create Base reaction and then fill it with its components.
    bounds, metabolites = _reaction_information(string=reaction)
    new_reaction = Reaction(id=reaction_id, name=reaction_name)
    for identifier, coefficient in metabolites.items():
        # Either get from model, or retrieve it.
        try:
            metabolite = model.metabolites.get_by_id(identifier)
        except KeyError:
            # _get_metabolite will also search for the metabolite under a
            # different name.
            compartment = identifier[-1:]
            new_identifier = identifier[:-2]
            try:
                data_dict = get_data(
                    directory=directory,
                    identifier=new_identifier,
                    database=database,
                    debug_level=10,
                )
                metabolite = _get_metabolite(
                    metabolite_dict=data_dict,
                    compartment=compartment,
                    model=model,
                )
            # It is necessary to build the metabolite.
            except (WrongParserError, HTTPError):
                debug_log.debug(
                    f'Metabolite "{identifier}" was identified as a manually '
                    + "curated reaction."
                )
                metabolite = _build_metabolite(
                    identifier=identifier,
                    formula="X",
                    name=identifier,
                    charge=0,
                    compartment=compartment,
                )
                msg = (
                    f'Manually-curated metabolite "{metabolite.id}" does not '
                    "have a chemical formula and charge. Please modify the "
                    'values manually. Charge set to 0, and formula to "X".'
                )
                debug_log.warning(msg=msg)
                warn(message=msg, category=UserWarning)

        new_reaction.add_metabolites({metabolite: coefficient})
        debug_log.debug(
            f'Metabolite "{metabolite.id}" added to Reaction '
            f'"{new_reaction.id}".'
        )
    # Set new bounds and check imbalance
    new_reaction.bounds = bounds
    check_imbalance(
        reaction=new_reaction,
        show_imbalance=show_imbalance,
        stop_imbalance=stop_imbalance,
    )
    return new_reaction


def _convert_string_reaction(
    line: str,
    model: Model,
    directory: Path,
    database: str,
    stop_imbalance: bool,
    show_imbalance: bool,
    replacement: dict = {},
    genome: str = None,
) -> Reaction:
    """
    Returns a Reaction from string. It can be either custom, or the identifier
    for a reaction in a database. The function will search for the reactions
    and its corresponding metabolites under other names inside the model and
    return it, if necessary.

    Syntax:

    :code:`reaction_identifier, compartment`

    For custom reactions:

    :code:`reaction_identifier, reaction_name | coefficient metabolite <->
    coefficient metabolite

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`4 OXYGEN-MOLECULE_c`**

    Else, include a model to retrieve identifiers from it.

    Reversibility will depend on the type of arrow that. Options are

    **`<--, -->, <=>, <->`**

    Args:
        line (str): String with custom reaction or identifier of reaction
        model (Model): Model to search for object if necessary.
        directory (Path): Path to directory, where data is located.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.
        genome (str): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie.
    """
    try:
        # Create custom reaction
        new_reaction = _reaction_from_string(
            line_string=line,
            directory=directory,
            database=database,
            model=model,
            stop_imbalance=stop_imbalance,
            show_imbalance=show_imbalance,
        )
    # If delimiter is not found, then it must be a reaction
    except NoDelimiter:
        # add reaction from root. Get only left part
        segment = (part.strip().rstrip() for part in line.split(","))
        identifier = next(segment)
        with suppress(KeyError):
            identifier = replacement[identifier]
        new_reaction = _obtain_reaction(
            model=model,
            identifier=identifier,
            directory=directory,
            database=database,
            compartment=next(segment),
            replacement=replacement,
            genome=genome,
            stop_imbalance=stop_imbalance,
            show_imbalance=show_imbalance,
        )
    return new_reaction


def _get_file_reactions(
    model: Model,
    filename: Path,
    stop_imbalance: bool,
    show_imbalance: bool,
    **kwargs,
) -> List[Reaction]:
    """
    Returns list with reactions from file. All reactions can be either created
    manually or retrieved from a database. For each reactions, its always
    checks for mass balance.

    Custom reactions:

    :code:`reaction_identifier, reaction_name | coefficient metabolite <->
    coefficient metabolite

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`4 OXYGEN-MOLECULE_c`**

    Else, include a model to retrieve identifiers from it. Reversibility will
    depend on the type of arrow that. Options are:

    **`<--, -->, <=>, <->`**

    From database:

    :code:`original_identifier, compartment`


    Args:
        model (Model): model to test
        filename (Path): location of the file with reaction information
        stop_imbalance (bool): If unbalanced reaction is found, stop process.
        show_imbalance (bool): If unbalanced reaction is found, show output.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie.

    Raises:
        FileNotFoundError: if file does not exists
    """
    # TODO: add mass balance check
    if not filename.exists():
        raise FileNotFoundError(
            f'File "{filename.name}" does not exist. '
            + "Check if the given path is correct."
        )
    with open(filename, "r") as f:
        new_reactions = list()
        for line in _read_lines(f=f):
            new_reactions.append(
                _convert_string_reaction(
                    line=line,
                    model=model,
                    stop_imbalance=False,
                    show_imbalance=True,
                    **kwargs,
                )
            )
    return new_reactions


def _ident_pathway(data_dict: dict) -> Generator:
    """
    Checks if given dictionary belongs to a pathway object. Returns its as
    generator if True. Else, raises WrongDataError.
    """
    if data_dict["TYPE"] == "Pathway":
        debug_log.debug(
            f"Object '{data_dict['ENTRY']}' identified as a pathway."
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
        show_imbalance (bool): If unbalanced reaction is found, show output.

    Returns:
       Generator: Reaction from dictionary.

    Raises:
        WrongDataError: If data does not include reaction information.
    """
    try:
        debug_log.debug(
            f"Object '{data_dict['ENTRY']}' identified as a reaction."
        )
        # First, create object and then yield it
        reaction = _get_reaction(
            data_dict=data_dict,
            directory=directory,
            replacement=replacement,
            database=database,
            compartment=compartment,
            model=model,
            stop_imbalance=stop_imbalance,
            show_imbalance=show_imbalance,
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
        # First, create object and then yield it
        metabolite = _get_metabolite(
            metabolite_dict=data_dict, compartment=compartment, model=model
        )
        debug_log.debug(f"Object '{metabolite.id}' identified as a metabolite")
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
    model_id: str = "universal",
    genome: str = None,
) -> Union[Reaction, Metabolite, dict]:
    """
    Creates and returns a COBRApy object based on given identifier and
    database. Identifier names will be formatted.

    .. hint:: Hyphens will become underscores. Double hyphens become single\
    underscores.

    Args:
        identifier (str): Original identifier of the database.
        directory (Path): Path to directory where data is stored.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names.
        compartment (str): Location of the object. In case of reaction, all
            metabolites will be included in the same location.
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}. Does not apply to
            pathways.

    Arguments for reactions:
        stop_imbalance (bool, optional): If an unbalanced reaction is found,
            stop the process. Defaults to False.
        show_imbalance (bool, optional): If an unbalanced reaction is found,
            print output. Defaults to True.
        model (Model, optional): Model to add search for translated metabolites
            or reactions. Defaults to an empty model.


    Special arguments for databases:
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Pathways are not available.
            Defaults to: "universal"
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            species involved. Genes will be obtained for this species.
            List available at https://www.genome.jp/kegg/catalog/org_list.html

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
        genome=genome,
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
    raise WrongDataError(
        f'Data of "{identifier}" cannot be identified properly. Please contact'
        + " the maintainers of CobraMod"
    )


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
        msg = (
            f'Metabolite "{member.id}" is already present in the model. '
            "Skipping addition."
        )
        debug_log.warning(msg=msg)
        warn(message=msg, category=UserWarning)


def add_metabolites(model: Model, obj: Any, database=None, **kwargs):
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
        model (Model): Model to be expanded and searched for metabolites.
        obj: A Path; a list with either strings or Metabolite objects,
            or a single string. See syntax above.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names. Defaults
            to None (This is useful for custom metabolites).

    Keyword Arguments:
        directory (Path): Path to directory where the data is located.

    Raises:
        WrongSyntax (from str): If the syntax is not followed correctly as
            mentioned above.
        ValueError: If Keyword Arguments are missing.
        FileNotFoundError (from Path): if file does not exists
    """
    try:
        # In case of a Path
        if isinstance(obj, Path):
            # These variable will raise KeyError if no kwargs are passed.
            directory = kwargs["directory"]
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
                "Given object is not the type mentioned in the docstrings. "
                + "Please use a string, a list or a pathlib 'Path' object."
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
        msg = (
            f'Reaction "{member.id}" is already present in the model. '
            "Skipping addition."
        )
        debug_log.warning(msg=msg)
        warn(message=msg, category=UserWarning)


def add_reactions(
    model: Model,
    obj: Any,
    database=None,
    stop_imbalance: bool = False,
    show_imbalance: bool = True,
    **kwargs,
):
    """Adds given object into the model. The options are:

     - Path: A file with components. E. g. :
        Path.cwd().joinpath("file_with_names.txt")
     - List[Reactions]: A list with regular Reactions
     - str: Either the identifier with its corresponding compartment or a
     string with all components. This applies for the Path option. E.g. :

        :code:`reaction_identifier, compartment`

        For custom reactions

        :code:`reaction_identifier, reaction_name | coefficient metabolite <->
        coefficient metabolite

        Identifiers of metabolites have to end with an underscore and a
        compartment:

            E.g. **`4 OXYGEN-MOLECULE_c`**

     - List[str]: A list with multiple str with the mentioned syntax.

    Args:
        model (Model): Model to expand and search for reactions.
        obj: A Path; a list with either strings or Reaction objects,
            or a single string. See syntax above.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names. Defaults
            to None (Useful for custom reactions).
        stop_imbalance (bool): If an unbalanced reaction is found,
            stop the process. Defaults to False.
        show_imbalance (bool): If an unbalanced reaction is found, print
            output. Defaults to True.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        replacement (dict): Original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.
        genome (str): Exclusive for KEGG. Abbreviation for the species
            involved. Genes will be obtained for this species.

    Raises:
        WrongSyntax (from str): If the syntax is not followed correctly as
            mentioned above.
        ValueError: If Keyword Arguments are missing.
        FileNotFoundError (from Path): If the file does not exists.

    """
    try:
        # In case of a Path
        if isinstance(obj, Path):
            # These variable will raise KeyError if no kwargs are passed.
            directory = kwargs["directory"]
            # Defaults
            replacement = kwargs.pop("replacement", dict())
            genome = kwargs.pop("genome", None)
            new_reactions = _get_file_reactions(
                model=model,
                filename=obj,
                directory=directory,
                database=database,
                replacement=replacement,
                stop_imbalance=stop_imbalance,
                show_imbalance=show_imbalance,
                genome=genome,
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
            # Defaults
            replacement = kwargs.pop("replacement", dict())
            genome = kwargs.pop("genome", None)
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
                    stop_imbalance=stop_imbalance,
                    show_imbalance=show_imbalance,
                    genome=genome,
                )
                for line in obj
            ]
        # Raise error if wrong
        else:
            raise WrongDataError(
                "Given object is not the type mentioned in the docstrings. "
                + "Please use a string, a list or a pathlib 'Path' object."
            )
        # Otherwise, it must be a list with Metabolites.
        __add_reactions_check(model=model, reactions=new_reactions)
    except KeyError:
        raise ValueError("Keyword Arguments are missing for given object.")
