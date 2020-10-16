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
        # CASE 1: simple retrieval from biocyc
        test_dict = get_data(
            directory=dir_data, database="META", identifier="WATER"
        )
        self.assertIsInstance(obj=test_dict, cls=dict)
        # CASE 2: simple retrieval from AraCyc
        # TODO: Case 2
        # CASE 3: Simple retrieval from KEGG
        test_dict = get_data(
            directory=dir_data, database="KEGG", identifier="C00001"
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
