from cobramod.debug import debug_log
from typing import Any, Dict


def _p_compound(root: Any) -> dict:
    """
    Parses given xml root into a dictionary for metabolite from biocyc
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    try:
        formula = root.find(
            "./*/cml/*/formula").attrib["concise"].replace(" ", "")
        charge = int(root.find("./*/cml/molecule").attrib["formalCharge"])
    # Must be Protein
    except AttributeError:
        formula = "X"
        charge = 0
        debug_log.warning(f'Biocyc ID "{identifier}" treated as Protein')
    # For names only
    try:
        name = root.find("./*/cml/*").attrib["title"]
    except AttributeError:
        name = identifier
    # Temporal fixup
    if formula == "":
        formula = "X"
    temp_dict = {
        "TYPE": root.find("*[@frameid]").tag,
        "ENTRY": identifier,
        "NAME": name,
        "FORMULA": formula,
        "CHARGE": charge
    }
    return temp_dict


def _p_metabolites(root: Any) -> dict:
    """
    Parses the data of given root and returns a dictionary with the
    information.
    """
    meta_dict: Dict[str, float] = dict()
    left_metabolites = root.findall("./Reaction/left")
    right_metabolites = root.findall("./Reaction/right")
    for meta in left_metabolites:
        try:
            coef = float(meta.find("coefficient").text) * -1
        except AttributeError:
            coef = -1  # default
        try:
            meta_identifier = meta.find(
                "*/[@frameid]").attrib["frameid"].strip().rstrip()
        except AttributeError:
            raise AttributeError('Reaction cannot find participants')
        meta_dict[meta_identifier] = coef
    for meta in right_metabolites:
        try:
            coef = float(meta.find("coefficient").text)
        except AttributeError:
            coef = 1  # default
        try:
            meta_identifier = meta.find(
                "*/[@frameid]").attrib["frameid"].strip().rstrip()
        except AttributeError:
            raise AttributeError('Reaction cannot find participants')
        meta_dict[meta_identifier] = coef
    return meta_dict


def _check_direction(root: Any) -> tuple:
    """
    Veifies that the direction of the reactions is the same as stated in
    the root file.
    """
    # Reversible <->
    text = root.find("*/reaction-direction").text
    if "REVERSIBLE" in text:
        bounds = (-1000, 1000)
    elif "RIGHT-TO-LEFT" in text:
        bounds = (-1000, 0)
    elif "LEFT-TO-RIGHT" in text:
        bounds = (0, 1000)
    else:
        raise Warning(
           f'Reversibility for '
           f'"{root.find("*[@frameid]").attrib["frameid"]}" could not '
           f'be found')
    return bounds


def _p_reaction(root: Any) -> dict:
    """
    Parses the root file and return the data in a dictionary
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    try:
        name = root.find(
            "*[@ID]/enzymatic-reaction/*/common-name").text
    except AttributeError:
        name = identifier
    temp_dict = {
        "TYPE": root.find("*[@frameid]").tag,
        "ENTRY": identifier,
        "NAME": name,
        "EQUATION": _p_metabolites(root=root),
        "BOUNDS": _check_direction(root=root)
    }
    return temp_dict


def _get_unsorted_pathway(root: Any) -> tuple:
    """
    Returns a dictionary with sequences (edges) for a graph; and a set
    with all participants (vertex).
    """
    reaction_dict = dict()
    reaction_set = set()
    for rxn_line in root.findall("*/reaction-ordering"):
        current = rxn_line.find("*/[@frameid]").attrib["frameid"]
        prior = rxn_line.find("predecessor-reactions/").attrib["frameid"]
        # If the direction of keys and values changes, then
        # many reactions would get lost. This way, even with
        # multiple compounds, pathways remain
        # NOTE: check for behaviour
        # Replacing values produces cuts with are needed to avoid cyclic
        # No reactions are missed
        reaction_dict[current] = prior
        # TODO: add information
        reaction_set.add(current)
        reaction_set.add(prior)
    name = root.find("Pathway").attrib["frameid"]
    # If dictinary has only one element, raise an error.
    if len(reaction_dict) == 0:
        raise NotImplementedError(
            'Path has only a reaction. Add separately')
    debug_log.debug(
        f'Dictionary for pathway "{name}" succesfully created')
    return reaction_dict, reaction_set


def _p_pathway(root: Any) -> dict:
    """
    Return dictionary with data for given xml root.
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    reaction_dict, reaction_set = _get_unsorted_pathway(root=root)
    temp_dict = {
        "TYPE": root.find("*[@frameid]").tag,
        "NAME": root.find("*[@frameid]").attrib["frameid"],
        "ENTRY": identifier,
        "PATHWAY": reaction_dict,
        "SET": reaction_set
        }
    return temp_dict
