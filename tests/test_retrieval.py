#!/usr/bin/env python3
"""
Unit-test for data retrieval and Data version configuration. It checks that
the files are loaded and saved properly
"""

import unittest
from pathlib import Path

from cobra import __version__ as cobra_version

import cobramod.retrieval as cmod_retrieval
from cobramod import __version__ as cmod_version
from cobramod.debug import change_to_debug
from cobramod.parsing.db_version import DataVersionConfigurator

change_to_debug()

dir_data = Path(__file__).resolve().parent.joinpath("data")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class RetrievalTesting(unittest.TestCase):
    @classmethod
    def setUp(cls):
        data_conf = DataVersionConfigurator()
        data_conf.ignore_db_versions = True

    def test_get_data(self):
        dir_data.joinpath("META", "WATER.xml").unlink(True)
        self.assertFalse(dir_data.joinpath("META", "WATER.xml").exists())

        # Simple retrieval from Biocyc (metabolite)
        test_data = cmod_retrieval.get_data(
            directory=dir_data, database="META", identifier="WATER"
        )
        self.assertEqual(test_data.path, dir_data.joinpath("META", "WATER.xml"))
        self.assertTrue(test_data.path.exists())

        # simple retrieval from Biocyc (reaction)
        test_data = cmod_retrieval.get_data(
            directory=dir_data, database="ARA", identifier="GLYOXII-RXN"
        )
        self.assertEqual(
            test_data.path, dir_data.joinpath("ARA", "GLYOXII-RXN.xml")
        )
        self.assertTrue(test_data.path.exists())
        self.assertTrue(
            dir_data.joinpath("ARA", "GENES", "GLYOXII-RXN_genes.xml").exists(),
        )

        # Simple retrieval from KEGG (metabolite)
        test_data = cmod_retrieval.get_data(
            directory=dir_data, database="KEGG", identifier="C00001"
        )
        self.assertTrue(test_data.path.exists())
        self.assertTrue(dir_data.joinpath("KEGG", "C00001.txt").exists())

        # Simple retrieval from KEGG (Reaction)
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            database="KEGG",
            identifier="R02736",
            genome="eco",
        )
        self.assertEqual(
            test_data.path, dir_data.joinpath("KEGG", "R02736.txt")
        )
        self.assertTrue(test_data.path.exists())
        self.assertTrue(
            dir_data.joinpath("KEGG", "GENES", "R02736_genes.txt").exists(),
        )

        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="R08618", database="KEGG"
        )
        self.assertEqual(
            test_data.path, dir_data.joinpath("KEGG", "R08618.txt")
        )
        self.assertTrue(test_data.path.exists())
        self.assertTrue(
            dir_data.joinpath("KEGG", "GENES", "R08618_genes.txt").exists(),
        )

        # Simple retrieval from KEGG (Pathway)
        test_data = cmod_retrieval.get_data(
            directory=dir_data, database="KEGG", identifier="M00001"
        )
        self.assertEqual(
            test_data.path, dir_data.joinpath("KEGG", "M00001.txt")
        )
        self.assertTrue(test_data.path.exists())

        # Simple retrieval from BIGG (metabolite)
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            database="BIGG",
            identifier="accoa_c",
            model_id="e_coli_core",
        )
        self.assertEqual(
            test_data.path,
            dir_data.joinpath("BIGG", "e_coli_core", "accoa_c.json"),
        )
        self.assertTrue(test_data.path.exists())

        # Simple retrieval from BIGG (reaction)
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            database="BIGG",
            identifier="CS",
            model_id="e_coli_core",
        )
        self.assertEqual(
            test_data.path, dir_data.joinpath("BIGG", "e_coli_core", "CS.json")
        )
        self.assertTrue(test_data.path.exists())

        # Retrieval from PMN
        test_data = cmod_retrieval.get_data(
            directory=dir_data, database="pmn:CORN", identifier="RXN-11501"
        )
        self.assertEqual(
            test_data.path, dir_data.joinpath("PMN", "CORN", "RXN-11501.xml")
        )
        self.assertTrue(test_data.path.exists())
        self.assertTrue(
            dir_data.joinpath(
                "PMN", "CORN", "GENES", "RXN-11501_genes.xml"
            ).exists(),
        )


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
