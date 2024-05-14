import filecmp
import tempfile
import unittest
from pathlib import Path

import cobra

from cobramod.test import textbook_biocyc
from cobramod.visualization.escher import EscherIntegration


class TestEscherIntegration(unittest.TestCase):
    def test_create(self):
        escher = EscherIntegration()
        self.assertIsInstance(escher, EscherIntegration)

    def test_save_html(self):
        escher = EscherIntegration()

        escher.map_json = (
            Path(__file__).parent.joinpath("escher_ref_map.json").read_text()
        )
        escher_html_expected = Path(__file__).parent.joinpath(
            "escher_html_ref.html"
        )

        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory) / "escher.html"
            escher.save_html(file)

            self.assertTrue(file.exists())
            self.assertTrue(filecmp.cmp(file, escher_html_expected))
