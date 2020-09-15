#!/usr/bin/env python3
from pathlib import Path
import xml.etree.ElementTree as ET
from cobra import Metabolite, Model, Reaction
from typing import Union, TextIO, Iterator
import requests
import logging

# Creating corresponding Logs
# Format
debug_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# Handler
debug_handler = logging.FileHandler("debug.log", mode="a+")
debug_handler.setFormatter(debug_formatter)
# Log
debug_log = logging.getLogger("debug_log")
debug_log.setLevel(logging.DEBUG)
# GenLog.ad
debug_log.addHandler(debug_handler)


def _define_base_dir(directory: Path, database: str) -> Path:
    """
    Returns Path object for given database. If directory does not exist.
    It will be created.

    Args:
        directory (Path): Parent directory.
        database (str): Name of database. Options: "META", "ARA".

    Returns:
        Path: Path object for database.
    """
    if directory.joinpath(database).exists():
        return directory.joinpath(database)
    else:
        directory.joinpath(database).mkdir()
        return directory.joinpath(database)


def get_xml_from_biocyc(
        directory: Path, identifier: str, database: str,
        **kwargs) -> ET.Element:
    """
    Searchs in given parent directory if data is located in their respective
    database directory. If not, data will be retrievied from the corresponding
    database. Returns root of given identifier.

    Args:
        directory (Path): Path to directory where data is located.
        identifier (str): identifier for given database.
        database (str): Name of database. Options: "META", "ARA".

    Raises:
        Warning: If object is not available in given database
        FileNotFoundError: If parent directory is not found

    Returns:
        ET.Element: root of XML file
    """
    if directory.exists():
        data_dir = _define_base_dir(
            directory=directory, database=database)
        filename = data_dir.joinpath(f'{identifier}.xml')
        debug_log.debug(f'Searching "{identifier}" in directory "{database}"')
        try:
            root = ET.parse(str(filename)).getroot()
            debug_log.debug('Found')
            return root
        except FileNotFoundError:
            debug_log.warning(
                f'"{identifier}" not found in directory "{database}".')
            # Retrieve from URL
            url_text = (
                f'https://websvc.biocyc.org/getxml?{database}:{identifier}')
            debug_log.debug(f'Searching in {url_text}')
            r = requests.get(url_text)
            if r.status_code == 404:
                msg = f'"{identifier}" not available in "{database}"'
                debug_log.error(msg)
                raise Warning(msg)
            else:
                root = ET.fromstring(r.text)  # defining root
                tree = ET.ElementTree(root)
                tree.write(str(filename))
                debug_log.debug(
                    f'Object found and saved in directory "{database}".')
                return root
    else:
        msg = "Directory not found"
        debug_log.error(msg)
        raise FileNotFoundError(msg)


def _create_meta_from_string(line_string: str) -> Metabolite:
    """
    Creates a Metabolite object based on a string.
    The string must follow the syntax:
    'formatted identifier, name , compartment, chemical_formula,
    molecular_charge'

    Args:
        line_string (str): string with information

    Raises:
        TypeError: if no string is identifier
        IndexError: if format is invalid

    Returns:
        Metabolite: New Metabolite with given information.
    """
    if not isinstance(line_string, str):
        raise TypeError('Argument must be a str')
    line = [part.strip().rstrip() for part in line_string.split(",")]
    try:
        meta_id = line[0]
        meta_name = line[1]
        meta_comp = line[2]
        meta_formula = line[3]
        meta_charge = line[4]
    except IndexError:
        raise IndexError(
            'Given line is invalid. Format is: id, name, compartment, '
            'chemical_formula, molecular_charge')
    return Metabolite(
        id=meta_id, name=meta_name, compartment=meta_comp, charge=meta_charge,
        formula=meta_formula)


def create_meta_from_root(
    root: Union[ET.Element, str], compartment: str = "c",
        **kwargs) -> Metabolite:
    """
    Creates a Metabolite object base on a root from a XML file.
    If not database is found, it will automatically search in the database
    "META".

    Args:
        root (Union[ET.Element, str]): root of XML file or identifier for
            specific database
        compartment (str, optional): [description]. Defaults to "c".

    Raises:
        TypeError: if given root is invalid

    Returns:
        Metabolite: New Metabolite based on root
    """
    if isinstance(root, str):
        try:
            root = get_xml_from_biocyc(identifier=root, **kwargs)
        except Warning:
            kwargs["database"] = "META"
            root = get_xml_from_biocyc(identifier=root, **kwargs)
    if not isinstance(root, ET.Element):
        raise TypeError('Given root is not valid')
    id_base = root.find("*/[@frameid]").attrib["frameid"]
    try:
        formula = root.find(  # chemical formula
            "./*/cml/*/formula").attrib["concise"].replace(" ", "")
        # obtaining molecular charge
        charge = int(
            root.find("./*/cml/molecule").attrib["formalCharge"])
    # TODO: distinguish oxidized and reduced members
    except AttributeError:  # must be an enzyme
        formula = "X"
        charge = 0
    id_base = id_base.replace("-", "_") + "_" + compartment
    try:
        name = root.find("*/cml/*").attrib["title"]
    except AttributeError:
        name = id_base
    return Metabolite(
        id=id_base, formula=formula, name=name, charge=charge,
        compartment=compartment)


def _read_lines(f: TextIO) -> Iterator:
    """
    Reads Text I/O and returns iterator of line that are not comments nor
    blanks spaces
    """
    for line in f:
        line = line.strip()
        if line.startswith("#"):
            continue
        if not line:  # blank
            continue
        yield line


def _has_root_name(line: str) -> bool:
    """
    Returns whether given line includes a unformatted identifier.
    """
    # TODO test some names
    line_separated = [part.strip().rstrip() for part in line.split(",")]
    return "_" not in line_separated[0][-3:]


def _check_if_meta_in_model(metabolite: str, model: Model, **kwargs) -> bool:
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
            f'Metabolite "{metabolite.id}" was found in given model. Skipping')
    else:
        model.add_metabolites([metabolite])
        debug_log.info(f'Metabolite "{metabolite.id}" was added to model')


def _has_comma_separator(line: str) -> bool:
    """Returns True if given line has a comma"""
    return "," in line


def _get_name_compartment_string(line: str) -> list:
    """Returns list of words separated previously by a comma"""
    parts = [part.strip().rstrip() for part in line.split(",")]
    return parts


def add_meta_from_string(
        line: str, model: Model, replacement_dict: dict = {},
        **kwargs) -> Metabolite:
    """
    Transform a string into a Metabolite object and appends it into model.
    The Metabolite can be either custom or from a database. Returns new
    Metabolite object.

    For custom metabolite the syntax follows:
    'formatted identifier, name , compartment, chemical_formula,
    molecular_charge'

    For metabolite from root, syntax is:
    'identifier, compartment'

    Args:
        line (str): string with either custom metabolite information or
            metabolite identifier from specific database
        model (Model): model
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        comparment (str): location of metabolites. Defaults to cytosol "c"
    Returns:
        Metabolite: new created Metabolite object
    """
    if _has_root_name(line=line):
        # FIX: to gain perfomance, search for name and then create Metabolite
        if _has_comma_separator(line=line):
            parts = _get_name_compartment_string(line)
            root = parts[0]
            kwargs["compartment"] = parts[1]
        else:
            # This is specific for metabolites coming from reactions
            root = line
        try:
            newMeta = create_meta_from_root(
                root=replacement_dict[root], **kwargs)
            debug_log.warning(
                f'Metabolite "{root}" replaced by '
                f'"{replacement_dict[root]}".')
        except KeyError:
            newMeta = create_meta_from_root(
                root=root, **kwargs)
        _add_if_not_found_model(model=model, metabolite=newMeta)
    else:
        newMeta = _create_meta_from_string(line_string=line)
        _add_if_not_found_model(model=model, metabolite=newMeta)
    return newMeta


def add_meta_from_file(model: Model, filename: Path, **kwargs):
    """
    Creates new Metabolites specified in given file. Syntax is mentioned in
    function 'add_meta_from_string'

    Args:
        model (Model): model to test
        filename (Path): location of the file with metabolites
        **kwargs: same as 'add_meta_from_string'

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
        lines = list(_read_lines(f=f))
        [add_meta_from_string(
            line=single, model=model, **kwargs) for single in lines]


def _create_base_reaction(
        root: ET.Element, compartment: str = "c", **kwargs) -> Reaction:
    """
    From given root from a XML file, creates and returns a base reaction.
    Location, identifier and name are formatted for the Reaction object.

    Args:
        root (Union[ET.Element, str]): root of XML file or identifier for
            specific database
        comparment (str, Optional): location of the reactions to take place.
            Defaults to cytosol "c"

    Keyword Arguments:
        comparment (str): location of the reactions to take place. Defaults to
            cytosol "c"
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.


    Returns:
        Reaction: Reaction object with basic information
    """
    base_id = root.find(
        "*/[@frameid]").attrib["frameid"].replace("--", "-").replace("-", "_")
    base_id = f'{base_id}_{compartment}'
    try:
        base_name = root.find(
            "*[@ID]/enzymatic-reaction/*/common-name").text
    except AttributeError:
        base_name = base_id
    return Reaction(id=base_id, name=base_name)


def _create_sides_for_reaction(
        model: Model, root: ET.Element, temp_reaction: Reaction,
        side: str = "left", **kwargs) -> Reaction:
    """
    For given root, creates either the reactant or product side of the reaction
    and adds it to given temporal Reaction object. Metabolites are retrieved
    from model or from given database.

    Args:
        model (Model): model to test
        root (ET.Element): root of XML file or identifier for
            specific database
        temp_reaction (Reaction): Base Reaction object
        side (str, optional): [description]. Defaults to "left".
        **kwargs: Same as '_create_base_reaction'
    Raises:
        TypeError: if side option not a string
        Warning: if side option is invalid
        AttributeError: if participants cannot be found.

    Returns:
        Reaction: Reaction object with new product or reactant side
    """
    if side == "left":
        MULTIPLIER = -1
        side_metabolites = root.findall("./Reaction/left")
    elif side == "right":
        MULTIPLIER = 1
        side_metabolites = root.findall("./Reaction/right")
    elif not isinstance(side, str):
        raise TypeError('Argument side is not valid')
    else:
        raise Warning('Only options are "right" and "left"')
    for meta in side_metabolites:
        try:
            coef = int(meta.find("coefficient").text) * MULTIPLIER
        except AttributeError:
            coef = MULTIPLIER  # default
        # NOTE: add replacement Dictionary
        try:
            meta_id = meta.find(
                "*/[@frameid]").attrib["frameid"].strip().rstrip()
        except AttributeError:
            raise AttributeError('Reaction cannot find participants')

        temp_metabolite = add_meta_from_string(
            line=meta_id, model=model, **kwargs)
        temp_reaction.add_metabolites({
            model.metabolites.get_by_id(temp_metabolite.id): coef})
    return temp_reaction


def _check_change_direction_reaction(reaction: Reaction, root: ET.Element):
    """
    Verifies that the direction of the reactions is the same as stated in the
    root file.
    """
    # Reversible <->
    text = root.find("*/reaction-direction").text
    if "REVERSIBLE" in text:
        reaction.bounds = (-1000, 1000)
    elif "RIGHT-TO-LEFT" in text:
        reaction.bounds = (-1000, 0)
    elif "LEFT-TO-RIGHT" in text:
        reaction.bounds = (0, 1000)


def build_reaction_from_xml(
        root: Union[ET.Element, str], **kwargs) -> Reaction:
    """
    Creates a Reactions Object from given root. Metabolites are searched in
    given model, otherwise retrieved from a specified database. If metabolite
    is not found, it will be search in "META"

    Args:
        root (Union[ET.Element, str]): root of XML file or identifier for
            specific database

    Keyword Arguments:
        model (Model): model to look up for metabolites
        comparment (str): location of the reactions to take place. Defaults to
            cytosol "c"
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.

    Returns:
        Reaction: New Reaction object
    """
    if isinstance(root, str):
        try:
            root = get_xml_from_biocyc(identifier=root, **kwargs)
        except Warning:
            kwargs["database"] = "META"
            root = get_xml_from_biocyc(identifier=root, **kwargs)
    # retrieve name and create simple reaction
    new_reaction = _create_base_reaction(root=root, **kwargs)
    # add left side
    new_reaction = _create_sides_for_reaction(
        root=root, temp_reaction=new_reaction, side="left", **kwargs)
    # direction
    _check_change_direction_reaction(reaction=new_reaction, root=root)
    # add right side
    new_reaction = _create_sides_for_reaction(
        root=root, temp_reaction=new_reaction, side="right", **kwargs)
    return new_reaction


def _check_if_reaction_in_model(reaction_id, model: Model) -> bool:
    """
    Returns whether reaction is found in model
    """
    return reaction_id in [reaction.id for reaction in model.reactions]


def _add_if_no_reaction_model(model: Model, reaction: Reaction):
    """
    Adds given Reaction objecto into model if this was not in the model,
    """
    if _check_if_reaction_in_model(model=model, reaction_id=reaction.id):
        debug_log.warning(
            f'Reaction "{reaction.id}" was found in given model. Skipping')
    else:
        model.add_reactions([reaction])
        debug_log.info(f'Reaction "{reaction.id}" was added to model')


def add_reaction_from_root(
        model: Model, root: Union[ET.Element, str],
        replacement_dict: dict = {}, **kwargs):
    """
    Creates a Reaction object from given object and adds it to given model if
    not found in the model

    Args:
        model (Model): model to test
        root (Union[ET.Element, str]): root of XML file or identifier for
            specific database
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.

    Keyword Arguments:
        comparment (str): location of the reactions to take place. Defaults to
            cytosol "c"
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".

    Raises:
        TypeError: If model is invalid
        TypeError: If given root is invalid
    """
    # validating variables
    if not isinstance(model, Model):
        raise TypeError('Model is not valid')
    if isinstance(root, str):
        try:
            root = get_xml_from_biocyc(
                identifier=replacement_dict[root], **kwargs)
            debug_log.warning(
                f'Replacing "{root}" with "{replacement_dict[root]}"')
        except KeyError:
            root = get_xml_from_biocyc(identifier=root, **kwargs)
    if not isinstance(root, ET.Element):
        raise TypeError('Given root is not valid')
    new_reaction = build_reaction_from_xml(
        root=root, model=model, **kwargs)
    _add_if_no_reaction_model(model=model, reaction=new_reaction)


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

    Syntax follows:
    'id_metabolite1: coefficient, id_metabolite2:coefficient ...
    id_metaboliteX: coefficient'

    Identifier has to end with an underscore and a compartment:
    E.g OXYGEN-MOLECULE_c: -1

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
        raise TypeError('Line format is wrong')
    tmpMetaDict = dict()
    for single in string_list:
        single = [x.strip().rstrip() for x in single.split(":")]
        meta_name = single[0]
        try:
            meta_value = float(single[1])
            tmpMetaDict[meta_name] = meta_value
        except ValueError:
            raise ValueError(f'Coefficient might be missing for {meta_name}')
    return tmpMetaDict


def create_custom_reaction(line_string: str, **kwargs) -> Reaction:
    """
    For given string, which includes the information of the Reaction and its
    metabolites. If metabolites are not in given model, it will be retrieved
    from specified database.

    The Syntax should follow:
    'reaction_identifier, reaction_name | metabolite_identifier1: coefficient,
    metabolite_identifier2:coefficient,...,metabolite_identifierX: coefficient'

    Identifier has to end with an underscore and a compartment:
    E.g OXYGEN-MOLECULE_c: -1

    Args:
        line_string (str): string with information

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.

    Raises:
        IndexError: if not identfier '|' is not found
        Warning: if identifier has a wrong format

    Returns:
        Reaction: new custom Reaction object
    """
    line_string = [x.strip().rstrip() for x in line_string.split("|")]
    rxn_id_name = [x.strip().rstrip() for x in line_string[0].split(",")]
    try:  # no blanks
        string_list = [x.strip().rstrip() for x in line_string[1].split(",")]
        metaDict = _build_dict_for_metabolites(string_list=string_list)
    except IndexError:  # wrong format
        raise IndexError(
            f'No delimiter "|" found for {rxn_id_name[0].split(",")[0]}.')
    if len(rxn_id_name) > 1 and rxn_id_name != [""]:
        rxn_id = rxn_id_name[0]
        rxnName = rxn_id_name[1]
    elif rxn_id_name == [""]:  # blank space
        raise Warning('Format is wrong. ID not detected')
    elif len(rxn_id_name) == 1:
        rxn_id = rxnName = rxn_id_name[0]
    new_reaction = Reaction(id=rxn_id, name=rxnName)
    for meta, coef in metaDict.items():
        if _check_if_meta_in_model(metabolite=meta, **kwargs):
            new_reaction.add_metabolites({
                kwargs["model"].metabolites.get_by_id(meta): coef})
        else:
            # FIX: avoid double creation of metabolites
            add_meta_from_string(
                line=meta[:-2], compartment=meta[-1],
                **kwargs)
            new_reaction.add_metabolites({
                kwargs["model"].metabolites.get_by_id(
                    meta.replace("-", "_")): coef})
    # FIXME: making all reversible
    new_reaction.bounds = (-1000, 1000)
    return new_reaction


def _add_reaction_line_to_model(line: str, model: Model, **kwargs):
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
            line_string=line, model=model, **kwargs)
        _add_if_no_reaction_model(model=model, reaction=new_reaction)
    else:
        # add reaction from root. Get only left part
        # !! not recognizing compartments
        line = [part.strip().rstrip() for part in line.split(",")]
        # !!
        try:
            add_reaction_from_root(
                model=model, root=line[0], compartment=line[1], **kwargs)
        except Warning:
            kwargs["database"] = "META"
            add_reaction_from_root(
                model=model, root=line[0], compartment=line[1], **kwargs)


def add_reactions_from_file(
        model: Model, filename: Path, **kwargs):
    """
    Adds new reactions to given Model. All reactions can be either created
    manually or retrieved from a database. For each reactions, its always
    checks for mass balance.

    The Syntax should follow:
    'reaction_identifier, reaction_name | metabolite_identifier1: coefficient,
    metabolite_identifier2:coefficient,...,metabolite_identifierX: coefficient'

    Identifier has to end with an underscore and a compartment:
    E.g OXYGEN-MOLECULE_c: -1

    Args:
        model (Model): model to test
        filename (Path): location of the file with reaction information

    Keyword Arguments:
        directory (Path): Path to directory where data is located.
        database (str): Name of database. Options: "META", "ARA".
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers.

    Raises:
        TypeError: If model is invalid
        FileNotFoundError: is file does not exists
    """
    # TODO: add mass balance check
    if not isinstance(model, Model):
        raise TypeError("Model given is not a valid")
    if not filename.exists():
        raise FileNotFoundError
    with open(filename, "r") as f:
        lines = list(_read_lines(f=f))
        [_add_reaction_line_to_model(
            line=single, model=model, **kwargs) for single in lines]


def check_mass_balance(
        model: Model, rxn_id: str, show_wrong: bool = True,
        stop_wrong: bool = False):
    """
    Verifies if given reaction is unbalanced in given model.

    Args:
        model (Model): model to test
        rxn_id (str): reaction identifier
        show_wrong (bool, optional): If unbalance is found, it will show the
            output. Defaults to True.
        stop_wrong (bool, optional): If unbalanace is found, raise a Warning.
            Defaults to False.

    Raises:
        Warning: if given reaction is unbalanced.
    """
    dict_balance = model.reactions.get_by_id(rxn_id).check_mass_balance()
    # Will stop if True
    if show_wrong is True and dict_balance != {}:
        msg = (
            f'Reaction unbalanced found at {rxn_id}. '
            f'Results to {dict_balance}. ')
        print(f'\n{msg}')
        if stop_wrong is True:
            msg += 'Stopping'
            raise Warning(msg)
        debug_log.warning(msg)
