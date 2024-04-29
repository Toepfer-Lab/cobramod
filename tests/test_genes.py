#!/usr/bin/env python3
"""Unittest for genes verification

In this module, the creation of genes for multiple function and their behavior
are checked
"""

import logging
import unittest
from pathlib import Path

import cobra.core as cobra_core
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.core.creation import create_object
from cobramod.debug import change_to_debug
from cobramod.parsing.db_version import DataVersionConfigurator

change_to_debug()
data_conf = DataVersionConfigurator()
data_conf.force_same_version = True

dir_data = Path(__file__).resolve().parent.joinpath("data")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class TestComplexGenes(unittest.TestCase):
    def test_KEGG_genome(self):
        # CASE: catch warning logs
        with self.assertLogs(level=logging.DEBUG) as cm:
            create_object(
                identifier="R02736",
                directory=dir_data,
                compartment="c",
                database="KEGG",
                genome="fake",
            )
            self.assertIn(
                "abbreviation as a specie. No genes will be added.",
                cm.output[0],
            )
        test_reaction = create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            genome="fake",
        )
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Given object is not a valid COBRApy Reaction")

        self.assertEqual(first=len(test_reaction.genes), second=0)

        # CASE: catch warning. No genome and thus, no genes
        with self.assertLogs(level=logging.DEBUG) as cm:
            create_object(
                identifier="R02736",
                directory=dir_data,
                compartment="c",
                database="KEGG",
            )
            self.assertIn(
                "Nothing was specified in argument ",
                cm.output[-1],
            )
        test_reaction = create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
        )
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Given object is not a valid COBRApy Reaction")

        self.assertEqual(first=len(test_reaction.genes), second=0)

        # CASE: regular case
        test_reaction = create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            genome="hsa",
        )
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Given object is not a valid COBRApy Reaction")

        self.assertEqual(first=len(test_reaction.genes), second=2)


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
