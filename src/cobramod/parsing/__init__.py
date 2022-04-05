#!/usr/bin/env python3
"""Sub-package parsing

This package handles the download of data and transformation into dictionaries
with can be later used for the creation of the corresponding COBRApy objects.
Classes included in the sub-package:

- BaseParser: Abstract class to download and parse data
- The corresponding child classes:
    - BiocycParser
    - PlantCycParser
    - KeggParser
    - BiggParser
    - SolCycParser

Read the corresponding documentation for each module and corresponding classes.
"""
from cobramod.parsing.base import BaseParser
from cobramod.parsing.biocyc import BiocycParser
from cobramod.parsing.plantcyc import PlantCycParser
from cobramod.parsing.bigg import BiggParser
from cobramod.parsing.kegg import KeggParser
from cobramod.parsing.solcyc import SolCycParser


__all__ = [
    "BaseParser",
    "BiocycParser",
    "PlantCycParser",
    "BiggParser",
    "KeggParser",
    "SolCycParser",
]
