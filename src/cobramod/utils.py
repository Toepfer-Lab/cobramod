"""CobraMod's helpful function

This module is intended to stored multiple functions, which cannot be grouped
to any module. In other words, these functions are general functions. An
example:

 - check_imbalance: Check for unbalanced reactions.
"""
from itertools import chain
from pathlib import Path
from re import match
from typing import Any, Generator, Iterable, Iterator, TextIO, Optional
from warnings import warn

import cobra.core as cobra_core
from cobra import DictList, Reaction

from cobramod.debug import debug_log
from cobramod.error import (
    PatternNotFound,
    UnbalancedReaction,
)

ARROWS: dict[str, tuple[int, int]] = {
    "<=>": (-1000, 1000),
    "<->": (-1000, 1000),
    "==>": (0, 1000),
    "-->": (0, 1000),
    "<==": (-1000, 0),
    "<--": (-1000, 0),
}


def build_metabolite(
    identifier: str, formula: str, name: str, charge: float, compartment: str
) -> cobra_core.Metabolite:
    """
    Returns a basic :class:`cobra.Metabolite`. Creation of the metabolite
    writes in the log a with a DEBUG level.
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
    metabolite = cobra_core.Metabolite(
        id=identifier,
        formula=formula,
        name=name,
        charge=charge,
        compartment=compartment,
    )
    debug_log.debug(f'Manually curated metabolite "{identifier}" was created.')
    return metabolite


def fix_name(name: str) -> str:
    """
    Returns new string, where hyphens are replaced to underscores
    """
    return name.replace("-", "_")


def check_imbalance(
    reaction: Reaction, stop_imbalance: bool, show_imbalance: bool
):
    """
    Verifies if given reaction is unbalanced in given model.

    Args:
        reaction (Reaction): Reaction object to examine.
        stop_imbalance (bool): If imbalance is found, stop process.
        show_imbalance (bool): If imbalance is found, show output.

    Raises:
        UnbalancedReaction: if given reaction is unbalanced.
    """
    try:
        dict_balance = reaction.check_mass_balance()
    except ValueError:
        dict()
    # Will stop if True
    if dict_balance != {}:
        msg = (
            f'Reaction "{reaction.id}" unbalanced. Following atoms are '
            + f"affected. Please verify:\n{dict_balance}"
        )
        if stop_imbalance:
            raise UnbalancedReaction(
                identifier=reaction.id, dict_balance=dict_balance
            )
        if show_imbalance:
            debug_log.warning(msg)
            warn(message=msg, category=UserWarning)


def get_key_dict(dictionary: dict, pattern: str) -> str:
    """
    From given pattern, return the first key of the dictionary that matches it
    or is a sub-string of it. It will raise a Warning if nothing found.
    """
    for key in dictionary.keys():
        if match(pattern=pattern, string=key):
            return str(key)
    raise PatternNotFound("No pattern was found for given dictionary.")


def read_lines(f: TextIO) -> Iterator:
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


def create_replacement(filename: Path) -> dict:
    """
    Creates a dictionary build from given file. Key are the first word until
    a colon ':' is found. The rest represents the value
    """
    replace_dict = dict()
    with open(file=str(filename), mode="r") as f:
        for line in read_lines(f=f):
            key, value = line.rsplit(sep=":")
            replace_dict[key.strip().rstrip()] = value.strip().rstrip()
    return replace_dict


def _replace_item(iterable: Iterable, replacement: dict = {}) -> Generator:
    """
    For an item in Iterable, replaces it for its corresponding value in
    given dictionary.

    Args:
        iterable (Iterable): sequence to modify
        replacement (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.

    Yields:
        Generator: Either original keys or the replacements
    """
    for item in iterable:
        if item in set(chain.from_iterable(replacement.keys())):
            debug_log.warning(f'Replacing "{item}" with "{replacement[item]}"')
            yield replacement[item]
        else:
            yield item


def _remove_item(iterable: Iterable, avoid_list: Iterable = []) -> Generator:
    """
    Returns Generator of items that are not in the avoid list.
    """
    for item in iterable:
        if item in avoid_list:
            debug_log.warning(f'Avoiding root for "{item}"')
        else:
            yield item


def compare_type(first: Any, second: Any):
    """
    Returns True is both objects are the same type, else raise TypeError.
    """
    if type(first) == type(second):
        return True
    else:
        raise TypeError("Given objects are not the same type.")


def compare_DictList(first: DictList, second: DictList) -> Generator:
    """
    Yields item identifiers from the first DictList that are not in the second
    DictList
    """
    compare_type(first=first[0], second=second[0])
    for item in first:
        if not second.has_id(id=item.id):
            yield item.id


def _save_diff(differences: dict) -> list:
    """
    Save differences in a list for later be used in a stdout or other kind of
    outputs.
    """
    output = list()
    if all([value == [] for value in differences.values()]):
        output.append("No differences!")
    else:
        output.append("Differences:")
        for key, value in differences.items():
            if value != []:
                output.append(f"{key.capitalize()}:")
                output += [f"- {identifier}" for identifier in value]
    return output


def write_to_file(sequences: Iterable, filename: Path):
    """
    Writes to given file all items from the iterable in separate lines.
    """
    with open(file=str(filename), mode="w+") as f:
        f.writelines(line + "\n" for line in sequences)


def _path_match(directory: Path, pattern: str) -> Path:
    """
    Returns first match as a Path object, given a pattern for a specific
    directory.
    """
    match = directory.glob(pattern=f"**/{pattern}.*")
    try:
        return next(match)
    except StopIteration:
        raise StopIteration(f"No file found with pattern '{pattern}'")


def find_intersection(
    dictlist: cobra_core.DictList, query: dict[str, str], revert: bool
) -> Optional[str]:
    """
    Return the first item from the intersection of a DictList and the values of
    a dictionary. The identifiers from the DictList can be reverted to their
    original. Method will raise a NoIntersectFound is no intersection is found.
    """
    # Format
    if revert:
        dict_set: set[str] = {
            item.id[:-2].replace("_", "-") for item in dictlist
        }
    # No format
    else:
        dict_set = {item.id for item in dictlist}
    # Error if nothing to pop, or empty Model
    try:
        common = dict_set.intersection(set(query.values()))
        return common.pop()
    except KeyError:
        return None


def is_compound(string: str) -> bool:
    try:
        float(string)

    except ValueError:
        if string[0] not in ["<", "-", "+"]:
            return True
    return False


def get_arrow_position(reaction_str: str) -> int:
    position = -1
    for item in ARROWS:
        position = reaction_str.find(item)

        if position != -1:
            break
    return position


def is_transport(reaction_str: str) -> bool:
    """
    Check if given reaction includes same participant in both sides,
    meaning that there is a tranport reaction. Prefix is not taken into
    consideration.
    """
    position = get_arrow_position(reaction_str)

    reactants = {
        item[2:].rstrip().strip()
        for item in reaction_str[:position].rstrip().strip().split(" ")
        if is_compound(item)
    }
    products = {
        item[2:].rstrip().strip()
        for item in reaction_str[position + 3 :].rstrip().strip().split(" ")
        if is_compound(item)
    }

    return len(reactants.intersection(products)) > 0


def convert_to_transport(reaction_str: str, compartment: str) -> str:
    """Modifies the reaction string such that the left side of the equation
    comes from the extracelullar compartment"""
    position = get_arrow_position(reaction_str)

    reactants = reaction_str[:position]
    products = reaction_str[position:]

    reactants = reactants.replace(f"_{compartment}", "_e")

    return reactants + products


def confirm_metabolite(
    model: cobra_core.Model, metabolites: list[cobra_core.Metabolite]
):
    """
    Checks if given metabolites are already in the model. If not, they will be
    added into the model.

    Args:
        model (Model): Model to extend.
        metabolites (List[Metabolites]): List with Metabolites.
    """
    # A Loop in necessary to log the skipped metabolites.
    for member in metabolites:
        found = model.metabolites.has_id(member.id)

        if not found:
            model.add_metabolites(metabolite_list=member)
            debug_log.info(f'Metabolite "{member.id}" was added to model.')
            continue

        msg = (
            f'Metabolite "{member.id}" is already present in the model. '
            "Skipping addition."
        )
        debug_log.warning(msg=msg)
        warn(message=msg, category=UserWarning)


def confirm_sink(model: cobra_core.Model, identifier: str):
    """
    Creates sink reactions for the metabolites that participate in given
    reaction. They can only be created if the corresponding metabolite only
    participates in one reaction. Otherwise, the metabolite can be somehow
    turnover
    """
    reaction = model.reactions.get_by_id(identifier)
    if not isinstance(reaction, cobra_core.Reaction):
        raise TypeError("Given object is not a COBRApy Reaction")

    metabolites: list[cobra_core.Metabolite] = (
        reaction.products + reaction.reactants
    )

    for metabolite in metabolites:
        amount = len(metabolite.reactions)

        if amount <= 1:
            # Warning is called in _inform_sink
            model.add_boundary(metabolite=metabolite, type="sink")


def add_reactions_to_model(
    model: cobra_core.Model, reactions: list[cobra_core.Reaction]
):
    """
    Check function that adds Reactions to given model if it does not
    contain the reaction. It logs the skipped reactions.
    """
    for member in reactions:
        if model.reactions.has_id(member.id):
            msg = (
                f'Reaction "{member.id}" is already present in the model. '
                "Skipping additions."
            )
            debug_log.warning(msg=msg)
            warn(message=msg, category=UserWarning)
            continue

        model.add_reactions([member])
        debug_log.info(f'Reaction "{member.id}" was added to model.')


def reaction_is_minimize(model: cobra_core.Model, identifier: str) -> bool:
    """
    Return whether given reaction should be minimized when optimizing
    """
    reaction = model.reactions.get_by_id(identifier)

    if not isinstance(reaction, cobra_core.Reaction):
        raise TypeError("Given object is not a COBRApy Reaction")

    lower = reaction.lower_bound
    upper = reaction.upper_bound

    if abs(lower) > upper:
        return True

    return False


def inform_new_sinks(model: cobra_core.Model, previous_sinks: set[str]):
    """
    Warns the user if sinks were created during the extension of the Model.
    i.e. Function compares a set with sink identifiers with the actual sinks
    in a model.
    """
    sinks: set[str] = {sink.id for sink in model.sinks if sink.id}.difference(
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
