"""
Check if models are available
"""

import unittest
from pathlib import Path

import pkg_resources


class ModuleTest(unittest.TestCase):
    def test_directory(self):
        self.assertTrue(
            Path(
                pkg_resources.resource_filename("cobra", "data/textbook.xml.gz")
            ).exists()
        )
        self.assertTrue(
            Path(
                pkg_resources.resource_filename(
                    "cobramod", "data/textbook_biocyc.sbml"
                )
            ).exists()
        )
        self.assertTrue(
            Path(
                pkg_resources.resource_filename(
                    "cobramod", "data/textbook_kegg.sbml"
                )
            ).exists()
        )


if __name__ == "__main__":
    unittest.main(verbosity=2, failfast=True)
