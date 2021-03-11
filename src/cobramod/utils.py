#!/usr/bin/env python3
"""CobraMod's helpful function

This module is intended to stored multiple functions, which cannot be grouped
to any module. In other words, these functions are general functions. An
example:

 - check_imbalance: Check for unbalanced reactions.
 - model_convert: Trasnform all Groups into proper Pathways.
"""
from itertools import chain
from pathlib import Path
from typing import TextIO, Iterator, Generator, Iterable, NamedTuple, Any, Dict
from re import match
from warnings import catch_warnings, simplefilter, warn

from cobra import Model, Reaction, DictList

from cobramod.debug import debug_log
from cobramod.error import (
    UnbalancedReaction,
    PatternNotFound,
    NoIntersectFound,
)
from cobramod.core.pathway import Pathway


class DataModel:
    """
    A class to store old values of metabolic models
    """

    def __init__(self, lists: Dict[str, list]):
        self.reactions = lists.get("reactions")
        self.metabolites = lists.get("metabolites")
        self.demands = lists.get("demands")
        self.exchanges = lists.get("exchanges")
        self.genes = lists.get("genes")
        self.groups = lists.get("groups")
        self.sinks = lists.get("sinks")

    @classmethod
    def from_model(cls, model: Model):
        data = {'reactions': model.reactions.list_attr("id"),
                'metabolites': model.metabolites.list_attr("id"),
                'demands': model.demands.list_attr("id"),
                'exchanges': model.exchanges.list_attr("id"),
                'genes': model.genes.list_attr("id"),
                'groups': model.groups.list_attr("id"),
                'sinks': model.sinks.list_attr("id")}
        return cls(data)

    def __sub__(self, other):
        data = {'reactions': [x for x in self.reactions if x not in set(other.reactions)],
                'metabolites': [x for x in self.metabolites if x not in set(other.metabolites)],
                'demands': [x for x in self.demands if x not in set(other.demands)],
                'exchanges': [x for x in self.exchanges if x not in set(other.exchanges)],
                'genes': [x for x in self.genes if x not in set(other.genes)],
                'groups': [x for x in self.groups if x not in set(other.groups)],
                'sinks': [x for x in self.sinks if x not in set(other.sinks)]}

        return DataModel(lists=data)

    def __str__(self):
        output = ['Reactions:\n',
                  str(self.reactions), "\n",
                  "Metabolites:\n",
                  str(self.metabolites), "\n",
                  "Exchange:\n",
                  str(self.exchanges), "\n",
                  "Demand:\n",
                  str(self.demands), "\n",
                  "Sinks:\n",
                  str(self.sinks), "\n",
                  "Genes:\n",
                  str(self.genes), "\n",
                  "Groups:\n",
                  str(self.groups), "\n"]

        return ''.join(output)


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
            f"Reaction '{reaction.id}' unbalanced. "
            f"Results to {dict_balance}. "
        )
        if stop_imbalance:
            raise UnbalancedReaction(reaction=reaction)
        if show_imbalance:
            debug_log.warning(msg)
            warn(message=msg, category=UserWarning)


def get_DataList(model: Model) -> DataModel:
    """
    Retrieve all DictList in given model and returns them as a
    :func:`cobramod.utils.DataModel` object.
    """
    # FIXME: copying a model is counterproductive. Attribute Group does not
    # have method copy(). Check what can be done
    copy_model = model.copy()
    dict_arguments = dict()
    with catch_warnings():
        # This is to avoid DeprecationWarning when using dir with Model
        simplefilter(action="ignore")
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


def _compare(model: Model, comparison: DataModel) -> dict:
    """
    Compares given model with the data from a DataModel object and returns
    a dictionary with differences.
    """
    difference = dict()
    with catch_warnings():
        # This is avoid DeprecationWarning when using dir with Model
        simplefilter(action="ignore")
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


def get_diff(model: Model, comparison: DataModel) -> list:
    """
    Gets the difference between give model and a
    :func:`cobramod.utils.DataModel` object and returns a list with
    differences.

    Args:
        model (Model): Model to obtain information
        comparison (DataModel): Object with data from a previous model. Check
            :func:`cobramod.utils.get_DataList` for usage.

    Returns:
        list: Differences between model and DataModel.

    """
    return _save_diff(differences=_compare(model=model, comparison=comparison))


def write_to_file(sequences: Iterable, filename: Path):
    """
    Writes to given file all items from the iterable in separate lines.
    """
    with open(file=str(filename), mode="w+") as f:
        f.writelines(line + "\n" for line in sequences)


def get_basic_info(model: Model) -> list:
    """
    Returns as a list the information of the model. The order of the items,
    represents the order for printing.
    """
    return [
        "Summary:",
        f"Model identifier: {model.id}",
        "Model name:",
        str(model.name),
        "Reactions:",
        str([reaction.id for reaction in model.reactions]),
        "Metabolites:",
        str([metabolite.id for metabolite in model.metabolites]),
        "Exchanges:",
        str([exchange.id for exchange in model.exchanges]),
        "Demands:",
        str([demand.id for demand in model.demands]),
        "Sinks:",
        str([sink.id for sink in model.sinks]),
        "Genes:",
        str([gene.id for gene in model.genes]),
        "Groups:",
        str([group.id for group in model.groups]),
    ]


def check_to_write(
        model: Model,
        summary: bool,
        filename: Path,
        old_values: DataModel,
):
    """
    Checks if summary should be saved in a filename.

    Args:
        model (Model): model with recent changes.
        summary (bool): If True, summary will be saved in given filename.
        filename (Path): Path with the location of the filename
        old_values (DataModel): Object with data from previous model. Use
            method :func:`cobramod.utils.get_DataList`.
        basic_info (list): Basic information of model. Use method
            :func:`cobramod.utils.get_basic_info`

    Raises:
        TypeError: if model is empty

    """
    try:
        if summary:

            new_values = DataModel.from_model(model)

            output = [
                "Summary:",
                f"Model identifier: {model.id}",
                "Model name:",
                str(model.name),
                str(new_values),
                "Added:",
                str(new_values - old_values),
                "Removed:",
                str(old_values - new_values)
            ]

            write_to_file(sequences=output, filename=filename)
        else:
            pass
    except TypeError:
        # To avoid empty models
        debug_log.error("Given model appears to be empty. Summary skipped")


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
    original. Method will raise a KeyError is no intersection is found.
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


def model_convert(model: Model):
    """
    Converts the all Group objects in given model to proper cobramod
    :func:`cobramod.pathway.Pathway`
    """
    for index, group in enumerate(iterable=model.groups):
        model.groups[index] = Pathway._transform(obj=group)
