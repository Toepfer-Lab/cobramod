from cobramod.debug import debug_log


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


class UnbalancedReaction(Exception):
    """
    Simple Error that should be raised if a reaction has wrong mass balance.
    """

    def __init__(self, reaction: str):
        """
        Args:
            reaction (str): identifier of the reaction.
        """
        debug_log.warning(f"Reaction {reaction} is not balanced.")


class NotInRangeError(Exception):
    """
    Simple Error that should be raised if, after adding a reaction, the
    optimized value, is not in range.
    """

    def __init__(self, reaction: str):
        """
        Args:
            reaction (str): identifier of the reaction.
        """
        debug_log.critical(
            f"Reaction '{reaction}' not in range. Check sinks manually."
        )
