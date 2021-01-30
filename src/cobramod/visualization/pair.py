from collections import UserDict
from contextlib import suppress

from cobramod.error import FoundInPairError


class PairDictionary(UserDict):
    """
    Dictionary that can include a pair (linked) dictionary. Both dictionaries
    cannot have the same keys. It is posible that only a dictionary is a
    PairDictionary, while the other a regular dictionary. The method
    :func:`cobramod.visualization.PairDictionary.set_pair' sets the pair
    for the dictionary.
    """

    def __init__(self, pair: UserDict = None, **kwargs):
        """
        Creates the regular dictionary and can set a pair dictionary if given.
        All keyword-arguments from a regular dictionary can be passed to this
        method.
        """
        self.pair = pair
        super().__init__(**kwargs)

    def __setitem__(self, key, value):
        """
        Set given value for given key in the PairDictionary. It will look
        in the keys of the pair and raise en FoundInPairError.
        """
        # Ignore None
        with suppress(AttributeError):
            if key in self.pair.keys():
                raise FoundInPairError(f"Pair dictionary has {key} as a key")
        super().__setitem__(key=key, item=value)

    def set_pair(self, pair: UserDict):
        """
        Set new pair for the dictionary. If will search for conflicts and raise
        a FoundInPairError if both PairDictionary have common keys.
        """
        # First check for problems and then set pair
        self._find_conflict(pair=pair)
        self.pair = pair

    def _find_conflict(self, pair: UserDict):
        """
        Takes the data of the PairDictionary and its pair and compares it. It
        will raise a FoundInPairError if both dictionaries share common keys.
        """
        first = set(self.data.keys())
        second = set(pair.keys())
        if first.intersection(second):
            raise FoundInPairError
