#!/usr/bin/env python3
from pathlib import Path
import xml.etree.ElementTree as ET
from cobra import Metabolite, Model, Reaction
from typing import Union, TextIO, Iterator
import requests
import logging
# Creating corresponding Logs
# Format
DebugFormatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# Handler
DebugHandler = logging.FileHandler("debug.log", mode="a+")
DebugHandler.setFormatter(DebugFormatter)
# Log
debug_log = logging.getLogger("debug_log")
debug_log.setLevel(logging.DEBUG)
# GenLog.ad
debug_log.addHandler(DebugHandler)


def _define_base_dir(directory: Path, database: str) -> Path:
    if directory.joinpath(database).exists():
        return directory.joinpath(database)
    else:
        directory.joinpath(database).mkdir()
        return directory.joinpath(database)


def get_xml_from_biocyc(
        directory: Path, identifier: str, database: str,
        **kwargs) -> ET.Element:
    # TODO: make differentation between database in files!! !!
    """Searchs in local DIR if given 'biocyc' .xml is found. If not, query
    Biocyc to retrive file. Returns root of given file.

    :param directory: Path to directory where data is located
    :type directory: Path
    :param identifier: Specific database identifier
    :type identifier: str
    :param database: Name for subdatabse. Options are: "META", "ARA"
    :type database: str, optional
    :raises Warning: If BioCyc Object ID is not found
    :raises FileNotFoundError: IF file is not located in directory
    :return: root of xml file
    :rtype: ET.Element
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
    """Reads Text I/O and returns iterator of line that are not comments nor
    blanks spaces

    :param f: Text input, e.g. a read file
    :type f: TextIO
    :yield: lines that are not comments or blanks
    :rtype: Iterator
    """
    for line in f:
        line = line.strip()
        if line.startswith("#"):
            continue
        if not line:  # blank
            continue
        yield line


def _has_root_name(line: str) -> bool:
    # TODO test some names
    line_separated = [part.strip().rstrip() for part in line.split(",")]
    return "_" not in line_separated[0][-3:]


def _check_if_meta_in_model(metaID, model: Model, **kwargs) -> bool:
    return metaID in [meta.id for meta in model.metabolites]


def _add_if_not_found_model(model: Model, metabolite: Metabolite):
    if _check_if_meta_in_model(model=model, metaID=metabolite.id):
        debug_log.warning(
            f'Metabolite "{metabolite.id}" was found in given model. Skipping')
    else:
        model.add_metabolites([metabolite])
        debug_log.info(f'Metabolite "{metabolite.id}" was added to model')


def _has_comma_separator(line: str) -> bool:
    """Returns True if given line has a comma"""
    return "," in line


def _get_name_compartment_string(line):
    """Returns list of words separated previously by a comma"""
    parts = [part.strip().rstrip() for part in line.split(",")]
    return parts


def add_meta_from_string(
        line: str, model: Model, replacement_dict: dict = {},
        **kwargs) -> Metabolite:
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
    """Creates new metabolites from given file and appends them into model.
    Metabolites can either be created with custom configuration or directly
    retrieved from BioCyc.
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
    return reaction_id in [reaction.id for reaction in model.reactions]


def _add_if_no_reaction_model(model: Model, reaction: Reaction):
    if _check_if_reaction_in_model(model=model, reaction_id=reaction.id):
        debug_log.warning(
            f'Reaction "{reaction.id}" was found in given model. Skipping')
    else:
        model.add_reactions([reaction])
        debug_log.info(f'Reaction "{reaction.id}" was added to model')


def add_reaction_from_root(
        model: Model, root: Union[ET.Element, str],
        replacement_dict: dict = {}, **kwargs):
    """Creates Reaction from given root. If no metabolites are found in
    Model, then rxn will search META
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
    """Returns true if given string has a vertical bar '|'

    :param line_string: line to check
    :type line_string: str
    :return: True if delimiter found
    :rtype: bool
    """
    return "|" in line_string


def _build_dict_for_metabolites(string_list: list) -> dict:
    """For given list with string, creates Dictionary where the keys are the IDs of
    metabolites and the values their corresponding values.
    Format should be:
    id_meta:value, id_meta2:value2...

    :param string_list: List with strings of "ID:coeffient"
    :type string_list: list
    :raises TypeError: if format is wrong
    :raises ValueError: If coeffient is missing (Wrong format)
    :return: Dictionary with new keys and values
    :rtype: dict
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
    """for given string which includes name of the Reaction and its components
    (metabolites), it creates a Reactions.

    The string should follow the syntax:
    rxn_id, rxnName | id_meta:value,id_meta2:value2....
    delimiter is a vertical bar '|'

    :param line_string: string with information
    :type line_string: str
    :raises IndexError: if no delimiter '|' is found (Wrong format)
    :raises Warning: if ID is not found (Wrong format)
    :return: new custom-created Reaction
    :rtype: Reaction
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
        if _check_if_meta_in_model(metaID=meta, **kwargs):
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
    """Adds new reactions to given Model. All reactions can be either created
    manually or retrieved from BioCyc. For each reactions, its always checks
    for mass balance. Unbalanced Reactions can either raise a Warning to
    be ignored.
    ### Format:
    A vertical bar "|" separates the reactions and metabolites in
    given file. It should follow the syntax:
    id_reactiom, name | id_meta_name:value, id_meta2_name:value, ...

    :param model: Model to add new reactions
    :type model: Model
    :param filename: Path of file with reactions
    :type filename: Path
    :param show_wrong: Print to console if reactions are unbalanced,
    defaults to True
    :type show_wrong: bool, optional
    :param stop_wrong: Raise Warning if reacion is unbalanced,
    defaults to False
    :type stop_wrong: bool, optional
    :raises Warning: If stopIfWrong is TRUE and reaction is unbalanced
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
    [summary]

    Args:
        model (Model): [description]
        rxn_id (str): [description]
        show_wrong (bool, optional): [description]. Defaults to True.
        stop_wrong (bool, optional): [description]. Defaults to False.

    Raises:
        Warning: [description]
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
