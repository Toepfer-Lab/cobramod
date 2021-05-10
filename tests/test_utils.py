#!/usr/bin/env python3
from contextlib import suppress
from logging import DEBUG
from pathlib import Path
import unittest

from cobra.core import Metabolite, Reaction

from cobramod.core.creation import create_object, add_reactions, get_data
from cobramod.debug import debug_log
from cobramod.error import NoIntersectFound
from cobramod.test import textbook_kegg
import cobramod.utils as ui

# Debug must be set in level DEBUG for the test
debug_log.setLevel(DEBUG)
# Setting directory for data
dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class UtilsTesting(unittest.TestCase):
    def test_check_imbalance(self):
        # Configuration
        test_reaction = create_object(
            directory=dir_data,
            identifier="RXN-11414",
            database="ARA",
            compartment="c",
            # In order to catch warning
            show_imbalance=False,
        )
        # CASE 1: Showing imbalance
        self.assertWarns(
            UserWarning,
            ui.check_imbalance,
            reaction=test_reaction,
            show_imbalance=True,
            stop_imbalance=False,
        )
        # CASE 2: Stopping at imbalance
        self.assertRaises(
            ui.UnbalancedReaction,
            ui.check_imbalance,
            reaction=test_reaction,
            show_imbalance=True,
            stop_imbalance=True,
        )

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

    def test_compare_type(self):
        # CASE 1: Native classes
        self.assertTrue(expr=ui.compare_type(first=dict(), second=dict()))
        self.assertRaises(TypeError, ui.compare_type, str(), dict())
        # CASE 2: COBRApy classes
        self.assertTrue(
            expr=ui.compare_type(first=Reaction(0), second=Reaction(1))
        )
        self.assertRaises(
            TypeError, ui.compare_type, Reaction(0), Metabolite(1)
        )

    def test_compare_DictList(self):
        # Preparation
        test_model = textbook_kegg.copy()
        add_reactions(
            model=test_model,
            directory=dir_data,
            database="KEGG",
            obj="R00894,c",
            replacement={},
        )
        # CASE 1: Different types
        test_generator = ui.compare_DictList(
            first=test_model.metabolites, second=test_model.reactions
        )
        self.assertRaises(TypeError, next, test_generator)
        # CASE 2: Regular Metabolites
        test_list = list(
            ui.compare_DictList(
                first=test_model.metabolites, second=textbook_kegg.metabolites
            )
        )
        self.assertEqual(first=len(test_list), second=3)
        # CASE 3: Regular Reactions
        test_list = list(
            ui.compare_DictList(
                first=test_model.reactions, second=textbook_kegg.reactions
            )
        )
        self.assertEqual(first=len(test_list), second=1)

    def test__print_differences(self):
        # CASE 1: regular dictionary
        test_dict = {"reactions": ["A", "B"], "metabolites": [1, 2]}
        test_diff = ui._save_diff(differences=test_dict)
        self.assertIn(member="- A", container=test_diff)
        # CASE 2: no differences
        test_dict = {"reactions": [], "metabolites": []}
        test_diff = ui._save_diff(differences=test_dict)
        self.assertIn(member="No differences!", container=test_diff)

    def test_save_to_file(self):
        # CASE 1: regular lines
        test_filename = dir_input.joinpath("summary.txt")
        with suppress(FileNotFoundError):
            test_filename.unlink()
        test_list = ["This should be first line", "This should second line"]
        ui.write_to_file(sequences=test_list, filename=test_filename)
        self.assertTrue(expr=test_filename.exists())
        with open(file=str(test_filename), mode="r") as e:
            self.assertEqual(first=2, second=sum(1 for line in e))
        test_filename.unlink()

    def test__path_match(self):
        # CASE 1a: regular match, KEGG compound
        test_path = ui._path_match(directory=dir_data, pattern="C00001")
        self.assertRegex(text=str(test_path), expected_regex="C00001.txt")
        # CASE 1b: regular match, KEGG reaction
        test_path = ui._path_match(directory=dir_data, pattern="R00114")
        self.assertRegex(text=str(test_path), expected_regex="R00114.txt")
        # CASE 2a: regular match, Biocyc compound
        test_path = ui._path_match(directory=dir_data, pattern="AMP")
        self.assertRegex(text=str(test_path), expected_regex="AMP.xml")
        # CASE 2b: regular match, Biocyc reaction
        get_data(
            directory=dir_data,
            identifier="GLUTAMINESYN-RXN",
            database="ARA",
            debug_level=10,
        )
        test_path = ui._path_match(
            directory=dir_data, pattern="GLUTAMINESYN-RXN"
        )
        self.assertRegex(
            text=str(test_path), expected_regex="GLUTAMINESYN-RXN.xml"
        )
        # CASE 3a: regular match, BIGG compound
        test_path = ui._path_match(directory=dir_data, pattern="accoa")
        self.assertRegex(text=str(test_path), expected_regex="accoa.json")
        # CASE 3b: regular match, BIGG compound in e_coli_core
        test_path = ui._path_match(directory=dir_data, pattern="accoa_c")
        self.assertRegex(text=str(test_path), expected_regex="accoa_c.json")

    def test__first_item(self):
        # CASE 1: Reactions Not found
        test_dict = {"A": "No", "B": "NotFound"}
        self.assertRaises(
            NoIntersectFound,
            ui._first_item,
            first=textbook_kegg.reactions,
            second=test_dict,
            revert=True,
        )
        # CASE 1: Metabolites
        test_dict = {"A": "C00058", "B": "NotFound"}
        test_string = ui._first_item(
            first=textbook_kegg.metabolites, second=test_dict, revert=True
        )
        self.assertEqual(first="C00058", second=test_string)
        # CASE 2: Reactions
        test_dict = {"A": "R00228", "B": "NotFound"}
        test_string = ui._first_item(
            first=textbook_kegg.reactions, second=test_dict, revert=True
        )
        self.assertEqual(first="R00228", second=test_string)


if __name__ == "__main__":
    unittest.main(verbosity=2)
