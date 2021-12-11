#!/usr/bin/env python3
"""CobraMod's helpful function

This module is intended to stored multiple functions, which cannot be grouped
to any module. In other words, these functions are general functions. An
example:

 - check_imbalance: Check for unbalanced reactions.
"""
from itertools import chain
from pathlib import Path
from re import match
from typing import TextIO, Iterator, Generator, Iterable, Any
from warnings import warn

from cobra import Reaction, DictList
from cobramod.debug import debug_log
from cobramod.error import (
    UnbalancedReaction,
    PatternNotFound,
    NoIntersectFound,
)


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
    dict_balance = reaction.check_mass_balance()
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


def create_replacement(filename: Path) -> dict:
    """
    Creates a dictionary build from given file. Key are the first word until
    a colon ':' is found. The rest represents the value
    """
    replace_dict = dict()
    with open(file=str(filename), mode="r") as f:
        for line in _read_lines(f=f):
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


def _first_item(first: DictList, second: dict, revert: bool) -> str:
    """
    Return the first item from the intersection of a DictList and the values of
    a dictionary. The identifiers from the DictList can be reverted to their
    original. Method will raise a NoIntersectFound is no intersection is found.
    """
    # Format
    if revert:
        dict_set = {item.id[:-2].replace("_", "-") for item in first}
    # No format
    else:
        dict_set = {item.id for item in first}
    # Error if nothing to pop, or empty Model
    try:
        common = dict_set.intersection(set(second.values()))
        return common.pop()
    except (KeyError, AttributeError):
        raise NoIntersectFound
