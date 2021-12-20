#!/usr/bin/env python3
"""Errors for CobraMod

This module creates special errors for CobraMod. Read each error for their
explanation.
"""
from cobramod.debug import debug_log


class UserInterruption(Exception):
    """
    Exception that is used if the input of a user interrupts the program flow.
    """

    pass


class NoDelimiter(Exception):
    """
    Simple Error that should be raised when a string does not include the
    delimiter "|"
    """

    def __init__(self, string: str):
        """
        Args:
            string (str): String with the format problem
        """
        msg = f'No delimiter "|" was found in \n{string}'
        super().__init__(msg)


class GraphKeyError(Exception):
    """
    Simple Error that should be raised when a value is missing as key in a
    graph
    """

    pass


class FoundInPairError(Exception):
    """
    Simple Error that is raised when the pair of a PairDictionary has already
    given key name.
    """

    pass


class NodeAttributeError(Exception):
    """
    Simple Error that should be raised when a
    :func:`cobramod.visualization.Node` cannot identify its type.
    """

    pass


class WrongParserError(Exception):
    """
    Simple Error that should be raised if a method cannot handle the parsing.
    """

    pass


class WrongDataError(Exception):
    """
    Simple Error that should be raised if COBRApy object is not the correct
    one.
    """

    pass


class WrongSyntax(Exception):
    """
    Simple Error that should be raised if a specific syntax is invalid. For
    instance, when strings do not have a certain number of items when splitted.
    """

    pass


class PatternNotFound(Exception):
    """
    Simple error that should be raised when a pattern is not found in a item.
    e.g. a substring.
    """

    pass


class NoIntersectFound(Exception):
    """
    Simple error that should be raised when no intersect is found between two
    sets.
    """

    pass


class NoGeneInformation(Exception):
    """
    Simple error that should be raised when given object has no gene
    information in the database.
    """

    pass


class AbbreviationWarning(Warning):
    """
    Simple Warning that should be raised when given abbreviation does not
    exists.
    """

    pass


class SuperpathwayWarning(Warning):
    """
    Simple Warning that should be raised if a pathway is identified as a
    Superpathway
    """

    pass


class UnbalancedReaction(Exception):
    """
    Simple Error that should be raised if a reaction has wrong mass balance.
    """

    def __init__(self, identifier: str, dict_balance: str):
        """
        Args:
            reaction (str): identifier of the reaction.
        """
        msg = (
            f'Reaction "{identifier}" unbalanced. Following atoms are '
            + f"affected. Please verify:\n{dict_balance}"
        )
        debug_log.critical(msg)
        super().__init__(msg)
