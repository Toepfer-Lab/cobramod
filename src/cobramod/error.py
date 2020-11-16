from cobramod.debug import debug_log


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
