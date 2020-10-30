#!/usr/bin/env python3
from pathlib import Path
from cobra import Metabolite, Model, Reaction
from typing import Union
import cobramod.mod_parser as par
from cobramod.utils import _read_lines
from cobramod.debug import debug_log
from contextlib import suppress
from collections import Counter


def _create_meta_from_string(line_string: str) -> Metabolite:
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


def _fix_name(name: str) -> str:
    """
    Replaces hyphens in given name to underscores. Double hyphens are
    transformed into single underscores
    """
    name = name.replace("--", "-")
    return name.replace("-", "_")


def build_metabolite(metabolite_dict: dict, compartment: str) -> Metabolite:
    """
    Builds and return a metabolite based on data saved as a dictionary.

    Args:
        metabolite_dict: dictionary with data of metabolite
        compartment: location of the metabolite

    Returns:
        Metabolite: New object based on a dictionary

    Raises:
        TypeError: If dictionary data does not have metabolite attributes
    """
    if not metabolite_dict["FORMULA"]:
        raise TypeError(
            "Given dictionary does not correspond to a Metabolite."
        )
    else:
        return Metabolite(
            id=f'{_fix_name(name=metabolite_dict["ENTRY"])}_{compartment}',
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


def meta_string_to_model(line: str, model: Model, **kwargs):
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
        database (str): Name of database. Options: "META", "ARA".

    Returns:
        Metabolite: New metabolite object

    Deprecated form: add_meta_from_string
    """
    # Kwargs belongs to query
    if line.count(",") > 1:
        # TODO: to gain perfomance, search for name and then create Metabolite
        new_metabolite = _create_meta_from_string(line_string=line)
    else:
        # Retrived from databse
        seqment = (part.strip().rstrip() for part in line.split(sep=","))
        # FIXME: include multiple databases
        # _identify_database()
        metabolite_dict = par.get_data(
            identifier=next(seqment), debug_level=10, **kwargs
        )
        new_metabolite = build_metabolite(
            metabolite_dict=metabolite_dict, compartment=next(seqment)
        )
    _add_if_not_found_model(model=model, metabolite=new_metabolite)
    return new_metabolite


def add_meta_from_file(model: Model, filename: Path, **kwargs):
    """
    Creates new Metabolites specified in given file. Syntax is mentioned in
    function :func:`cobramod.creation.meta_string_to_model`

    Args:
        model (Model): model to test
        filename (Path): location of the file with metabolites

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".

    Raises:
        TypeError: if model is invalid
        FileNotFoundError: if given file is not found
    """
    # checking validity of objects
    if not isinstance(model, Model):
        raise TypeError("Model given is not a valid")
    if not filename.exists():
        raise FileNotFoundError
    with open(filename, "r") as f:
        for line in _read_lines(f=f):
            try:
                meta_string_to_model(line=line, model=model, **kwargs)
            except Warning:
                # FIXME: add function to search for common missing databases
                kwargs["database"] = "META"
                meta_string_to_model(line=line, model=model, **kwargs)


def _return_duplicate(data_dict: dict) -> bool:
    """
    Check for the duplicate in a dictionary with prefixes-
    """
    sequence = Counter([item[2:] for item in data_dict.keys()])
    return sequence.most_common(1)[0][0]


def _build_reaction(
    data_dict: dict, compartment: str, replacement_dict: dict, **kwargs
) -> Reaction:
    """
    Creates Reactions object from given dictionary with data. Location of the
    reactions can be set with the argument 'compartment'.  Metabolites can be
    replaced by using the dictionary 'replacement_dict' with the following
    syntax:

    :code:`old_identifier: replacement`

    Args:
        data_dict (dict): dictionary with data of a Reaction.
        compartment (str): locations of the reactions
        replacement_dict (dict): original identifiers to be replaced.
            Values are the new identifiers.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".

    Returns:
        Reaction: New reaction based on dictionary

    Deprecated form: _create_base_reaction
    """
    reaction = Reaction(
        id=f'{_fix_name(name=data_dict["ENTRY"])}_{compartment}',
        name=data_dict["NAME"],
    )
    for identifier, coef in data_dict["EQUATION"].items():
        # Get rid of prefix r_ and l_
        identifier = identifier[2:]
        # TODO: add option to get metabolites from model
        # Only if found in replacemente
        with suppress(KeyError):
            identifier = replacement_dict[identifier]
            debug_log.warning(
                f'Metabolite "{identifier}" replaced by '
                f'"{replacement_dict[identifier]}".'
            )
        metabolite = par.get_data(
            identifier=identifier, debug_level=10, **kwargs
        )
        if (
            data_dict["TRANSPORT"]
            and coef < 0
            and _return_duplicate(data_dict=data_dict["EQUATION"])
            == identifier
        ):
            # FIX: temporary setting to extracellular
            metabolite = build_metabolite(
                metabolite_dict=metabolite, compartment="e"
            )
        else:
            metabolite = build_metabolite(
                metabolite_dict=metabolite, compartment=compartment
            )
        reaction.add_metabolites(metabolites_to_add={metabolite: coef})
        reaction.bounds = data_dict["BOUNDS"]
    return reaction


def _check_if_reaction_in_model(reaction_id, model: Model) -> bool:
    """
    Returns whether given reaction is found in model
    TODO: change to proper name
    """
    return reaction_id in [reaction.id for reaction in model.reactions]


def _add_if_no_reaction_model(model: Model, reaction: Reaction):
    """
    Adds given Reaction objecto into model if this was not in the model,
    """
    if _check_if_reaction_in_model(model=model, reaction_id=reaction.id):
        debug_log.warning(
            f'Reaction "{reaction.id}" was found in given model. Skipping'
        )
    else:
        model.add_reactions([reaction])
        debug_log.info(f'Reaction "{reaction.id}" was added to model')


def add_reaction(
    model: Model,
    identifier: str,
    directory: Path,
    database: str,
    compartment: str,
    replacement_dict: dict,
):
    """
    Creates and adds a Reaction object based on given identifier from given
    database. Metabolites can be replaced by using the dictionary
    replacement_dict with the following syntax:

    :code:`old_identifier: replacement`

    Args:
        model (Model): model
        identifier (str): identifier
        directory (Path): directory
        database (str): database
        compartment (str): compartment
        replacement_dict (dict): replacement_dict

    Deprecated form: add_reaction_from_root
    """
    data_dict = par.get_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=10,
    )
    reaction = _build_reaction(
        data_dict=data_dict,
        compartment=compartment,
        directory=directory,
        database=database,
        replacement_dict=replacement_dict,
    )
    _add_if_no_reaction_model(model=model, reaction=reaction)


def _has_delimiter(line_string: str) -> bool:
    """
    Returns true if given string has a vertical bar '|'
    """
    return "|" in line_string


def _build_dict_for_metabolites(string_list: list) -> dict:
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


def create_custom_reaction(
    line_string: str, directory: Path, database: str, model: Model = Model(0)
) -> Reaction:
    """
    For given string, which includes the information of the Reaction and its
    metabolites. If metabolites are not in given model, it will be retrieved
    from specified database.

    Syntax:

    :code:`reaction_identifier, reaction_name | metabolite_identifier1:
    coefficient,`
    :code:`metabolite_identifier2:coefficient, ..., metabolite_identifierX:
    coefficient`

    Identifier has to end with an underscore and a compartment:

    E.g **`OXYGEN-MOLECULE_c: -1`**

    Else, include a model to retrieve identifiers from it.

    Args:
        line_string (str): string with information

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.
        model (Model, Optional): A model to obtain metabolite Objects

    Raises:
        IndexError: if not identfier '|' is not found
        Warning: if identifier has a wrong format

    Returns:
        Reaction: New custom Reaction
    """
    segments = [x.strip().rstrip() for x in line_string.split("|")]
    rxn_id_name = [x.strip().rstrip() for x in segments[0].split(",")]
    try:  # no blanks
        string_list = [x.strip().rstrip() for x in segments[1].split(",")]
        meta_dict = _build_dict_for_metabolites(string_list=string_list)
    except IndexError:  # wrong format
        raise IndexError(
            f'No delimiter "|" found for {rxn_id_name[0].split(",")[0]}.'
        )
    if len(rxn_id_name) == 2:
        rxn_id, rxnName = rxn_id_name
    elif len(rxn_id_name) == 1 and rxn_id_name[0] != "":
        rxn_id = rxnName = rxn_id_name[0]
    else:  # blank space
        raise Warning(f"Wrong format for {segments}. ID not detected")
    new_reaction = Reaction(id=rxn_id, name=rxnName)
    for identifier, coef in meta_dict.items():
        try:
            # To avoid creation of metabolite if already in model
            metabolite = model.metabolites.get_by_id(identifier)
        except KeyError:
            compartment = identifier[-1:]
            identifier = identifier[:-2]
            data_dict = par.get_data(
                directory=directory,
                identifier=identifier,
                database=database,
                debug_level=10,
            )
            metabolite = build_metabolite(
                metabolite_dict=data_dict, compartment=compartment
            )
        new_reaction.add_metabolites({metabolite: coef})
    # FIXME: making all reversible
    new_reaction.bounds = (-1000, 1000)
    return new_reaction


def _add_reaction_line_to_model(
    line: str,
    model: Model,
    directory: Path,
    database: str,
    replacement_dict: dict = {},
    **kwargs,
):
    """
    From given string, it will identify if a custom Reaction or a Reaction
    from root can be created. It will build the reaction and adds it to
    given model

    Args:
        line (str): string with information
        model (Model): model to test
        **kwargs: same as in 'create_custom_reaction'
    """
    if _has_delimiter(line_string=line):
        # create custom reaction
        new_reaction = create_custom_reaction(
            line_string=line,
            directory=directory,
            database=database,
            model=model,
        )
        _add_if_no_reaction_model(model=model, reaction=new_reaction)
    else:
        # add reaction from root. Get only left part
        seqment = (part.strip().rstrip() for part in line.split(","))
        identifier = next(seqment)
        with suppress(KeyError):
            identifier = replacement_dict[identifier]
        # TODO: identify database
        add_reaction(
            model=model,
            identifier=identifier,
            directory=directory,
            database=database,
            compartment=next(seqment),
            replacement_dict=replacement_dict,
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
        replacement_dict (dict, optional): original identifiers to be replaced.
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


def check_mass_balance(
    model: Model,
    rxn_id: str,
    show_wrong: bool = True,
    stop_wrong: bool = False,
):
    """
    Verifies if given reaction is unbalanced in given model.

    Args:
        model (Model): model to test
        rxn_id (str): reaction identifier
        show_wrong (bool, optional): If unbalance is found, it wil/_buil show
            the output. Defaults to True.
        stop_wrong (bool, optional): If unbalanace is found, raise a Warning.
            Defaults to False.

    Raises:
        Warning: if given reaction is unbalanced.
    """
    dict_balance = model.reactions.get_by_id(rxn_id).check_mass_balance()
    # Will stop if True
    if show_wrong is True and dict_balance != {}:
        msg = (
            f"Reaction unbalanced found at '{rxn_id}'. "
            f"Results to {dict_balance}. "
        )
        print(f"\n{msg}")
        if stop_wrong is True:
            msg += "Stopping"
            raise Warning(msg)
        debug_log.warning(msg)


def create_object(
    identifier: str,
    directory: Path,
    database: str,
    compartment: str,
    replacement_dict: dict = {},
) -> Union[Reaction, Metabolite]:
    """
    Creates and returns COBRApy object based on given identifier and database.
    Identifier names will be formatted.

    .. hint:: Hyphens will become underscores. Double hyphens become single\
    underscores.

    Args:
        identifier (str): original identifier for database
        directory (Path): Path to directory where data is located.
        database (str): Name of the database. Options are: "META", "ARA",
            "KEGG"
        compartment (str): Location of the object. In case of reaction, all
            metabolites will be included in the same location
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.


    Returns:
        Union[Reaction, Metabolite]: Either Reaction on Metabolite object
    """
    data_dict = par.get_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=10,
    )
    # build_metabolite
    # FIX: Temporal solution. Pathways are missing
    # FIXME: logs information expressed twice
    try:
        debug_log.info(f"Metabolite for '{identifier}' identified")
        return build_metabolite(
            metabolite_dict=data_dict, compartment=compartment
        )
    except KeyError:
        debug_log.info(f"Reaction for '{identifier}' identified")
        return _build_reaction(
            data_dict=data_dict,
            directory=directory,
            replacement_dict=replacement_dict,
            database=database,
            compartment=compartment,
        )
