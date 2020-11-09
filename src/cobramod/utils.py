#!/usr/bin/env python3
from itertools import chain
from pathlib import Path
from typing import TextIO, Iterator, Generator, Iterable, NamedTuple, Any
from re import match

from cobra import Model, DictList

from cobramod.debug import debug_log


class WrongParserError(Exception):
    """
    Simple Error that should be raised if a method cannot handle the parsing
    """

    pass


class DataModel(NamedTuple):
    """
    A class to store old values of metabolic models
    """

    reactions: DictList
    metabolites: DictList
    demands: DictList
    exchanges: DictList
    genes: DictList
    groups: DictList
    sinks: DictList


def get_DataList(model: Model) -> DataModel:
    """
    Retrieve all DictList in given model and returns them as a
    :func:`cobramod.utils.DataModel` object.
    """
    # FIXME: copying a model is counterproductive. Attribute Group does not
    # have method copy()
    copy_model = model.copy()
    dict_arguments = dict()
    dict_list_names = [
        attribute
        for attribute in dir(model)
        if type(getattr(model, attribute)) == DictList
    ]
    for attribute in dict_list_names:
        item = getattr(copy_model, attribute)
        dict_arguments[attribute] = item
    return DataModel(**dict_arguments)


def get_key_dict(dictionary: dict, pattern: str) -> str:
    """
    From given pattern, return first key the dictionary that matches or is a
    sub-string.
    """
    for key in dictionary.keys():
        if match(pattern=pattern, string=key):
            return str(key)
    raise Warning("No pattern was found for given dictionary.")


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


def _replace_item(
    iterable: Iterable, replacement_dict: dict = {}
) -> Generator:
    """
    For an item in Iterable, replaces it for its corresponding value in
    given dictionary.

    Args:
        iterable (Iterable): sequence to modify
        replacement_dict (dict, optional): original identifiers to be replaced.
            Values are the new identifiers. Defaults to {}.

    Yields:
        Generator: Either original keys or the replacements
    """
    for item in iterable:
        if item in set(chain.from_iterable(replacement_dict.keys())):
            debug_log.warning(
                f'Replacing "{item}" with "{replacement_dict[item]}"'
            )
            yield replacement_dict[item]
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


def _compare(model: Model, comparison: DataModel) -> dict:
    """
    Compares given model with the data from a DataModel object  and returns
    a dictionary with differences.
    """
    difference = dict()
    dict_list_names = [
        attribute
        for attribute in dir(model)
        if type(getattr(model, attribute)) == DictList
    ]
    for dict_list in dict_list_names:
        item = getattr(model, dict_list)
        for key, value in comparison._asdict().items():
            try:
                if dict_list == key:
                    difference[dict_list] = list(
                        compare_DictList(first=item, second=value)
                    )
            except IndexError:
                difference[dict_list] = []
    return difference


def _print_differences(differences: dict):
    """
    Prints only differences
    """
    if all([item == [] for item in differences]):
        print("No differences!")
    else:
        print("Differences:")
        for key, value in differences.items():
            if value != []:
                print(f"{key.capitalize()}:")
                print(*[f"- {identifier}" for identifier in value], sep="\n")
