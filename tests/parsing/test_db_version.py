import logging
import shutil
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from pandas._testing import assert_frame_equal, assert_series_equal

from cobramod.parsing.db_version import DataVersionConfigurator

data_conf = DataVersionConfigurator()


class DataVersion(unittest.TestCase):
    directory: str

    @classmethod
    def setUp(cls):
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

    def test_multiple_directories(self):
        # CASE: Each data directory keeps its own versioning file
        with tempfile.TemporaryDirectory() as other_directory:
            data_conf.check_database_version(
                Path(other_directory), "pmn:META", "2.0"
            )
            data_conf.check_database_version(
                Path(self.directory), "kegg", "3.0"
            )

            other_versions = pd.read_csv(
                Path(other_directory) / "DatabaseVersions.csv", dtype=str
            )

        versions = pd.read_csv(
            Path(self.directory) / "DatabaseVersions.csv", dtype=str
        )

        assert_frame_equal(
            other_versions,
            pd.DataFrame({"orgid": ["pmn:META"], "version": ["2.0"]}),
        )
        assert_frame_equal(
            versions,
            pd.DataFrame(
                {"orgid": ["bigg", "kegg"], "version": ["1.0.0", "3.0"]}
            ),
        )

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
