#!/usr/bin/env python3
"""
Unit-test for data retrieval and Data version configuration. It checks that
the files are loaded and saved properly
"""

import logging
import shutil
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import cobramod.retrieval as cmod_retrieval
import pandas as pd
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.debug import change_to_debug
from cobramod.parsing.db_version import DataVersionConfigurator
from pandas._testing import assert_frame_equal, assert_series_equal

change_to_debug()
data_conf = DataVersionConfigurator()

dir_data = Path(__file__).resolve().parent.joinpath("data")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class DataVersion(unittest.TestCase):
    directory: str

    @classmethod
    def setUp(cls):
        # removing metadata belonging to other tests that affect this one
        data_conf.database_version = None

        cls.directory = tempfile.mkdtemp()
        data_conf.ignore_db_versions = False
        cls.versions = pd.DataFrame.from_dict(
            {"orgid": ["bigg"], "version": ["1.0.0"]}
        )

        cls.versions.to_csv(
            path_or_buf=str(cls.directory + "/DatabaseVersions.csv"),
            index=False,
        )

    @classmethod
    def tearDown(cls):
        shutil.rmtree(cls.directory)
        # Removing metadata
        data_conf.database_version = None

    def test_get_database_version(self):
        database = data_conf.get_database_version(self.directory)
        assert_frame_equal(database, self.versions)

    def test_get_local_databases(self):
        databases = data_conf.get_local_databases(self.directory)

        expected_series = pd.Series(["bigg"])
        expected_series.name = "orgid"
        assert_series_equal(databases, expected_series)

    def test_set_database_version(self):
        data_conf.set_database_version(Path(self.directory), "pmn:META", "1.0")

        database = data_conf.get_database_version(self.directory)

        self.versions = pd.concat(
            [
                self.versions,
                pd.DataFrame({"orgid": ["pmn:META"], "version": ["1.0"]}),
            ],
            ignore_index=True,
        )

        assert_frame_equal(database, self.versions)

    def test_check_database_version(self):
        (Path(self.directory) / "DatabaseVersions.csv").unlink()

        # new database
        data_conf.check_database_version(Path(self.directory), "bigg", "1.0.0")

        database = data_conf.get_database_version(self.directory)
        assert_frame_equal(database, self.versions)

        # database with correct version

        with warnings.catch_warnings(record=True) as warning_list:
            data_conf.check_database_version(
                Path(self.directory), "bigg", "1.0.0"
            )

        self.assertEqual(len(warning_list), 0)

        database = data_conf.get_database_version(self.directory)
        assert_frame_equal(database, self.versions)

        # BUG: if running discovery, this fails otherwise
        data_conf.force_same_version = False

        # database with incorrect version (with user input)
        with patch("builtins.input", return_value="y") as mock:
            with self.assertLogs(level=logging.DEBUG) as cm:
                data_conf.check_database_version(
                    Path(self.directory), "bigg", "2.0.0"
                )
                self.assertIn(
                    "do not match. Remote has version 2.0.0", cm.output[0]
                )

            args, _ = mock.call_args
            assert args == ("Ignore version mismatch? (Y)es (N)o: ",)

        database = data_conf.get_database_version(self.directory)
        assert_frame_equal(database, self.versions)

        # database with incorrect version (without user input)
        data_conf.ignore_db_versions = True
        with self.assertLogs(level=logging.DEBUG) as cm:
            data_conf.check_database_version(
                Path(self.directory), "bigg", "2.0.0"
            )
            self.assertIn(
                "do not match. Remote has version 2.0.0", cm.output[0]
            )

        database = data_conf.get_database_version(self.directory)
        assert_frame_equal(database, self.versions)


class RetrievalTesting(unittest.TestCase):
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
