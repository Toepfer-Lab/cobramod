#!/usr/bin/env python3
from logging import DEBUG
from pathlib import Path
import unittest

from cobramod.mod_parser import get_data
from cobramod.debug import debug_log

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
        # CASE 3b: Simple retrieval from KEGG (Reaction)
        test_dict = get_data(
            directory=dir_data, database="KEGG", identifier="R02736"
        )
        # CASE 3c: Simple retrieval from KEGG (Pathway)
        test_dict = get_data(
            directory=dir_data, database="KEGG", identifier="M00001"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Pathway")
        self.assertEqual(first=test_dict["ENTRY"], second="M00001")
        # CASE 4a: Simple retrieval from BIGG (metabolite)
        test_dict = get_data(
            directory=dir_data,
            database="BIGG",
            identifier="accoa_c",
            model_id="e_coli_core",
        )
        # CASE 4b: Simple retrieval from BIGG (reaction)
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        self.assertEqual(first=test_dict["ENTRY"], second="accoa")
        test_dict = get_data(
            directory=dir_data,
            database="BIGG",
            identifier="CS",
            model_id="e_coli_core",
        )
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        self.assertEqual(first=test_dict["ENTRY"], second="CS")

    def test_check_attribute(self):
        # CASE: Generic compounds:
        for test_dict in (
            get_data(
                directory=dir_data,
                database="META",
                identifier="Polyphosphates",
            ),
            get_data(directory=dir_data, database="KEGG", identifier="C00404"),
            get_data(directory=dir_data, database="KEGG", identifier="C00139"),
        ):
            self.assertEqual(first=test_dict["FORMULA"], second="X")


if __name__ == "__main__":
    unittest.main(verbosity=2)
