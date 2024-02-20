#!/usr/bin/env python3
import logging
import unittest
from contextlib import suppress
from pathlib import Path

import cobra.core as cobra_core
import cobramod.core.creation as cmod_core_creation
import cobramod.error as cmod_error
import cobramod.utils as ui
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.debug import change_to_debug
from cobramod.test import textbook_kegg

change_to_debug()

dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class UtilsTesting(unittest.TestCase):
    def test_check_imbalance(self):
        # Configuration
        test_reaction = cmod_core_creation.create_object(
            directory=dir_data,
            identifier="RXN-11414",
            database="ARA",
            compartment="c",
            # In order to catch warning
            show_imbalance=False,
        )
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")

        # CASE: Showing imbalance
        with self.assertLogs(level=logging.DEBUG) as cm:
            ui.check_imbalance(
                reaction=test_reaction,
                show_imbalance=True,
                stop_imbalance=False,
            )
            self.assertIn("unbalanced", cm.output[-1])

        # CASE: Stopping at imbalance
        self.assertRaises(
            cmod_error.UnbalancedReaction,
            ui.check_imbalance,
            reaction=test_reaction,
            show_imbalance=True,
            stop_imbalance=True,
        )

    def test_read_lines(self):
        # CASE 0: Comments and blank lines
        with open(file=dir_input.joinpath("test_reading_lines.txt")) as f:
            line = list(ui.read_lines(f=f))
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
            expr=ui.compare_type(
                first=cobra_core.Reaction("0"), second=cobra_core.Reaction("1")
            )
        )
        self.assertRaises(
            TypeError,
            ui.compare_type,
            cobra_core.Reaction("0"),
            cobra_core.Metabolite("1"),
        )

    def test_compare_DictList(self):
        # Preparation
        test_model = textbook_kegg.copy()
        cmod_core_creation.add_reactions(
            model=test_model,
            directory=dir_data,
            database="KEGG",
            obj=["R00894,c"],
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

    def test_save_to_file(self):
        # CASE 1: regular lines
        test_filename = dir_input.joinpath("summary.txt")
        with suppress(FileNotFoundError):
            test_filename.unlink()
        test_list = ["This should be first line", "This should second line"]
        ui.write_to_file(sequences=test_list, filename=test_filename)
        self.assertTrue(expr=test_filename.exists())
        with open(file=str(test_filename), mode="r") as e:
            self.assertEqual(first=2, second=sum(1 for _ in e))
        test_filename.unlink()

    def test_find_intersection(self):
        # CASE 1: Reactions Not found
        test_dict = {"A": "No", "B": "NotFound"}
        test_string = ui.find_intersection(
            dictlist=textbook_kegg.reactions,
            query=test_dict,
            revert=True,
        )
        self.assertIsNone(test_string)
        # CASE 1: Metabolites
        test_dict = {"A": "C00058", "B": "NotFound"}
        test_string = ui.find_intersection(
            dictlist=textbook_kegg.metabolites, query=test_dict, revert=True
        )
        self.assertEqual(first="C00058", second=test_string)
        # CASE 2: Reactions
        test_dict = {"A": "R00228", "B": "NotFound"}
        test_string = ui.find_intersection(
            dictlist=textbook_kegg.reactions, query=test_dict, revert=True
        )
        self.assertEqual(first="R00228", second=test_string)


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2)
