#!/usr/bin/env python3
import unittest
from pathlib import Path
from cobramod.debug import debug_log
from logging import DEBUG
import cobramod.utils as ui

debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")

if not dir_data.exists():
    dir_data.mkdir(parents=True)


class UtilsTesting(unittest.TestCase):
    def test__read_lines(self):
        # CASE 0: Comments and blank lines
        with open(file=dir_input.joinpath("test_reading_lines.txt")) as f:
            line = list(ui._read_lines(f=f))
        self.assertEqual(len(line), 4)

    def test_create_replacement(self):
        # CASE 1: regular dictionary with comments and empty spaces
        test_dict = ui.create_replacement(
            filename=dir_input.joinpath("test_create_replacement.txt")
        )
        self.assertIn(member="fum_e", container=test_dict.keys())
        self.assertIn(member="Glucopyranose_e", container=test_dict.values())


if __name__ == "__main__":
    unittest.main(verbosity=2)
