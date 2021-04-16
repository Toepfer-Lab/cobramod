#!/usr/bin/env python3
"""Unittest for genes verification

In this module, the creation of genes for multiple function and their behavior
are checked
"""
from logging import DEBUG
from pathlib import Path
from unittest import main, TestCase

from cobramod import create_object
from cobramod.debug import debug_log
from cobramod.error import AbbreviationWarning

# Debug must be set in level DEBUG for the test
debug_log.setLevel(DEBUG)
# Setting directory for data
dir_data = Path(__file__).resolve().parent.joinpath("data")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class TestComplexGenes(TestCase):
    def test_KEGG_genome(self):
        # CASE 1: catch FalseAbbreviation
        self.assertWarns(
            AbbreviationWarning,
            create_object,
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            genome="fake",
        )
        test_reaction = create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            genome="fake",
        )
        self.assertEqual(first=len(test_reaction.genes), second=0)
        # CASE 2: catch regular UserWarning. No genome and thus, no genes
        self.assertWarns(
            UserWarning,
            create_object,
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
        )
        test_reaction = create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
        )
        self.assertEqual(first=len(test_reaction.genes), second=0)
        # CASE 3: regular case
        test_reaction = create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            genome="hsa",
        )
        self.assertEqual(first=len(test_reaction.genes), second=2)


if __name__ == "__main__":
    main(verbosity=2)
