from typing import Optional

from attr import dataclass


@dataclass
class GenerellIdentifiers:
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    smiles: Optional[str] = None