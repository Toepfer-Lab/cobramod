#!/usr/bin/env python3
from logging import DEBUG
from pathlib import Path
import unittest

import cobramod.core.retrieval as md
from cobramod.debug import debug_log

# Debug must be set in level DEBUG for the test
debug_log.setLevel(DEBUG)
# Setting directory for data
dir_data = Path(__file__).resolve().parent.joinpath("data")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class RetrievalTesting(unittest.TestCase):
    def test_get_data(self):
        # CASE 1a: Simple retrieval from Biocyc (metabolite)
        test_dict = md.get_data(
            directory=dir_data, database="META", identifier="WATER"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 1b: simple retrieval from Biocyc (reaction)
        test_dict = md.get_data(
            directory=dir_data, database="ARA", identifier="GLYOXII-RXN"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        # CASE 3a: Simple retrieval from KEGG (metabolite)
        test_dict = md.get_data(
            directory=dir_data, database="KEGG", identifier="C00001"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 3b: Simple retrieval from KEGG (Reaction)
        test_dict = md.get_data(
            directory=dir_data,
            database="KEGG",
            identifier="R02736",
            genome="eco",
        )
        # CASE 3c: Simple retrieval from KEGG (Pathway)
        test_dict = md.get_data(
            directory=dir_data, database="KEGG", identifier="M00001"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Pathway")
        self.assertEqual(first=test_dict["ENTRY"], second="M00001")
        # CASE 4a: Simple retrieval from BIGG (metabolite)
        test_dict = md.get_data(
            directory=dir_data,
            database="BIGG",
            identifier="accoa_c",
            model_id="e_coli_core",
        )
        # CASE 4b: Simple retrieval from BIGG (reaction)
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        self.assertEqual(first=test_dict["ENTRY"], second="accoa")
        test_dict = md.get_data(
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
            md.get_data(
                directory=dir_data,
                database="META",
                identifier="Polyphosphates",
            ),
            md.get_data(
                directory=dir_data, database="KEGG", identifier="C00404"
            ),
            md.get_data(
                directory=dir_data, database="KEGG", identifier="C00139"
            ),
        ):
            self.assertEqual(first=test_dict["FORMULA"], second="X")

    def test_translate(self):
        # CASE 1: Regular compound KEGG
        md.get_data(directory=dir_data, database="KEGG", identifier="C00002")
        test_string = md.translate(
            directory=dir_data, target="C00002", database="CAS"
        )
        self.assertEqual(first=test_string, second="56-65-5")
        # CASE 2: Regular compound Biocyc
        md.get_data(directory=dir_data, database="ARA", identifier="AMP")
        test_string = md.translate(
            directory=dir_data, target="AMP", database="PUBCHEM"
        )
        self.assertEqual(first=test_string, second="15938965")
        # CASE 3: Regular compound BIGG
        test_string = md.translate(
            directory=dir_data, target="accoa_c", database="CHEBI"
        )
        self.assertEqual(first=test_string, second="CHEBI:13712")


if __name__ == "__main__":
    unittest.main(verbosity=2)
