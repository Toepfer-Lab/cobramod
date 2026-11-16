import filecmp
import json
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np

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

    def test_save_html_options(self):
        # CASE: Fluxes from COBRApy are numpy values and may contain None
        escher = EscherIntegration(
            map_json=Path(__file__)
            .parent.joinpath("escher_ref_map.json")
            .read_text(),
            reaction_data={"FUM": np.float64(5.0), "CS": None},
            reaction_scale=[
                {"type": "min", "color": "#0000ff", "size": 12},
                {
                    "type": "value",
                    "value": np.float64(0.0),
                    "color": "#ffffff",
                    "size": 16,
                },
            ],
            reaction_styles=["color", "size", "text"],
        )

        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory) / "escher.html"
            escher.save_html(file)
            html = file.read_text()

        # The options must be valid JSON, otherwise the map is not loaded
        match = re.search(r"\{\"reaction_styles\".*\}", html)
        self.assertIsNotNone(match)
        options = json.loads(match.group(0))

        self.assertEqual(options["reaction_data"], {"FUM": 5.0, "CS": None})
        self.assertEqual(
            options["reaction_scale"][1],
            {"type": "value", "value": 0.0, "color": "#ffffff", "size": 16},
        )
        self.assertIs(options["never_ask_before_quit"], False)
