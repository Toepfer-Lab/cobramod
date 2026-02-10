from types import NoneType

from pydantic import ConfigDict
from pydantic.dataclasses import dataclass

from typing import Dict, Union

@dataclass
class Unavailable:
    """
    Defines the absence of an ID from a database
    """

    def __str__(self):
        return "Unavailable"

    pass


class GenerellIdentifiers:
    """
    Identifiers can either be available eg. a string, unqueried defined as None or
     Unavailable if the Database did not return a results-
    """


    __inchi: Union[None, str, Unavailable] = None
    inchi_key: Union[None, str, Unavailable] = None
    smiles: Union[None, str, Unavailable] = None

    def __init__(self,
                 inchi: Union[None, str, Unavailable] = None,
                 inchi_key: Union[None, str, Unavailable] = None,
                 smiles: Union[None, str, Unavailable] = None,
                 ):
        super().__init__()

        self.inchi = inchi
        self.inchi_key = inchi_key
        self.smiles = smiles

    @property
    def inchi(self) -> Union[None, str, Unavailable]:
        return self.__inchi

    @inchi.setter
    def inchi(self, inchi: Union[None, str, Unavailable]):
        assert isinstance(inchi, (str, Unavailable, NoneType))

        if isinstance(inchi, str):
            assert inchi.startswith("InChI=")

        self.__inchi = inchi

    def weakEq(self, other):

        if self.inchi_key is None or other.inchi_key is None:
            pass
        elif self.inchi_key != other.inchi_key:
            return False


        if self.smiles is None or other.smiles is None:
            pass
        elif self.smiles != other.smiles:
            return False


        if self.inchi_key is None or other.inchi_key is None:
            pass
        elif self.inchi_key != other.inchi_key:
            return False

        return True

    def __iadd__(self, other):
        if self.inchi_key is None:
            self.inchi_key = other.inchi_key

        if self.smiles is None:
            self.smiles = other.smiles

        if self.inchi is None:
            self.inchi = other.inchi

        return self


    def empty(self):
        l = [self.smiles, self.inchi, self.inchi_key]
        if l.count(None) == len(l):
            return True
        else:
            return False

    @classmethod
    def fromAnnotation(cls, annotation: Dict[str, str]):
        smiles = annotation.get("smiles", None)
        inchi = annotation.get("inchi", None)
        inchi_key = annotation.get("inchi_key", None)

        return cls(
            smiles=smiles,
            inchi=inchi,
            inchi_key=inchi_key,
        )

    def to_dict(self):
        return {
            "smiles": str(self.smiles),
            "inchi": str(self.inchi),
            "inchi_key": str(self.inchi_key),
        }

    @classmethod
    def from_dict(cls, dict:Dict):

        smiles = dict["smiles"]
        inchi = dict["inchi"]
        inchi_key = dict["inchi_key"]

        if smiles == "Unavailable":
            smiles = Unavailable()
        elif smiles == "None":
            smiles = None

        if inchi_key == "Unavailable":
            inchi_key = Unavailable()
        elif inchi_key == "None":
            inchi_key = None

        if inchi == "Unavailable":
            inchi = Unavailable()
        elif inchi == "None":
            inchi = None

        return cls(
            smiles=smiles,
            inchi=inchi,
            inchi_key=inchi_key,
        )

    def anyNoneEntries(self) -> bool:

        if self.smiles is None:
            return True

        if self.inchi_key is None:
            return True

        if self.inchi is None:
            return True

        return False