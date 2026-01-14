from pydantic import ConfigDict
from pydantic.dataclasses import dataclass

from typing import Optional, Dict

from jedi.inference.gradual import annotation


@dataclass
class GenerellIdentifiers:
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    smiles: Optional[str] = None

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

    __pydantic_config__ = ConfigDict(validate_assignment=True)
