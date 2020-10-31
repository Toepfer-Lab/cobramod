#!/usr/bin/env python3
import unittest
from cobramod.mod_parser import get_data
from cobramod.debug import debug_log
from logging import DEBUG
from pathlib import Path

debug_log.setLevel(DEBUG)
dir_data = Path.cwd().joinpath("tests").joinpath("data")


class RetrievalTesting(unittest.TestCase):
    def test_get_data(self):
        # CASE 1a: Simple retrieval from Biocyc (metabolite)
        test_dict = get_data(
            directory=dir_data, database="META", identifier="WATER"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 1b: simple retrieval from Biocyc (reaction)
        test_dict = get_data(
            directory=dir_data, database="META", identifier="GLYOXII-RXN"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        # CASE 2: simple retrieval from AraCyc
        # TODO: Case 2
        # CASE 3a: Simple retrieval from KEGG (metabolite)
        test_dict = get_data(
            directory=dir_data, database="KEGG", identifier="C00001"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 3a: Simple retrieval from KEGG (Reaction)
        test_dict = get_data(
            directory=dir_data, database="KEGG", identifier="R02736"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        # CASE 4a: Simple retrieval from BIGG (metabolite)
        test_dict = get_data(
            directory=dir_data,
            database="BIGG",
            identifier="accoa_c",
            model_id="e_coli_core",
        )
        # CASE 4b: Simple retrieval from BIGG (reaction)
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        test_dict = get_data(
            directory=dir_data,
            database="BIGG",
            identifier="CS",
            model_id="e_coli_core",
        )
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")


if __name__ == "__main__":
    unittest.main(verbosity=2)
