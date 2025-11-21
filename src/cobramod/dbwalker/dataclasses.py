from dataclasses import dataclass
from typing import Optional


@dataclass
class GenerellIdentifiers:
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    smiles: Optional[str] = None
