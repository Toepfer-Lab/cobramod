"""
Check if models are available
"""

import unittest
from importlib import resources

import cobra

import cobramod


class ModuleTest(unittest.TestCase):
    def test_directory(self):
        self.assertTrue(
            resources.files(cobra).joinpath("data/textbook.xml.gz").is_file()
        )
        self.assertTrue(
            resources.files(cobramod)
            .joinpath("data/textbook_biocyc.sbml")
            .is_file()
        )
        self.assertTrue(
            resources.files(cobramod)
            .joinpath("data/textbook_kegg.sbml")
            .is_file()
        )


if __name__ == "__main__":
    unittest.main(verbosity=2, failfast=True)
