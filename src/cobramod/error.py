#!/usr/bin/env python3
"""Errors for CobraMod

This module creates special errors for CobraMod. Read each error for their
explanation.
"""
from cobramod.debug import debug_log

from cobra import Reaction, Metabolite


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


class NotInRangeError(Exception):
    """
    Simple Error that should be raised if, after adding a reaction, the
    optimized value, is not in range.
    """

    def __init__(self, reaction: Reaction):
        """
        Args:
            reaction (Reaction): COBRApy Reaction that cannot pass the
                non-zero flux test
        """
        # Count the number of reactions for metabolites and check which
        # metabolite has problem with its turnover.
        problem = []
        metabolite: Metabolite
        for metabolite in reaction.metabolites:
            reactions = [
                reaction.id
                for reaction in metabolite.reactions
                if "DM_" not in reaction.id or "SK_" not in reaction.id
            ]
            if len(reactions) == 1:
                problem.append(metabolite.id)
        msg = (
            f'The following reaction "{reaction.id}" failed the non-zero flux '
            "test multiple times. Flux values are below solver tolerance. "
        )
        if problem:
            msg += (
                "Please curate manually by adding reactions that enable "
                f'turnover of metabolites: {", ".join(problem)}'
            )
        else:
            msg += (
                "It it possible that one of the metabolites participates in "
                "a cycle. For example, NADP to NADPH and viceversa. Please "
                "make sure that those reactions have a correct equation and "
                "that their metabolites can be turnover."
            )
        debug_log.critical(msg)
        super().__init__(msg)
