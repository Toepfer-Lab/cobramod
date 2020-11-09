#!/usr/bin/env python3
from logging import DEBUG
from pathlib import Path
import unittest

from cobra import Metabolite, Reaction

from cobramod.debug import debug_log
from cobramod.test import textbook_kegg
from cobramod.creation import add_reaction
import cobramod.utils as ui

debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")

if not dir_data.exists():
    dir_data.mkdir(parents=True)


class UtilsTesting(unittest.TestCase):
    def test_get_DataList(self):
        # CASE 1: regular retrieval
        test_object = ui.get_DataList(model=textbook_kegg)
        self.assertIsInstance(obj=test_object, cls=ui.DataModel)

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
        add_reaction(
            model=test_model,
            compartment="c",
            directory=dir_input,
            database="KEGG",
            identifier="R00114",
            replacement_dict={},
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
        self.assertEqual(first=len(test_list), second=6)
        # CASE 3: Regular Reactions
        test_list = list(
            ui.compare_DictList(
                first=test_model.reactions, second=textbook_kegg.reactions
            )
        )
        self.assertEqual(first=len(test_list), second=1)

    def test__compare(self):
        test_model = textbook_kegg.copy()
        test_data = ui.get_DataList(model=test_model)
        add_reaction(
            model=test_model,
            compartment="c",
            directory=dir_input,
            database="KEGG",
            identifier="R00114",
            replacement_dict={},
        )
        test_dict = ui._compare(model=test_model, comparison=test_data)
        self.assertEqual(first=len(test_dict["metabolites"]), second=6)
        self.assertEqual(first=len(test_dict["reactions"]), second=1)

    def test__print_differences(self):
        test_model = textbook_kegg.copy()
        test_data = ui.get_DataList(model=test_model)
        add_reaction(
            model=test_model,
            compartment="c",
            directory=dir_input,
            database="KEGG",
            identifier="R00114",
            replacement_dict={},
        )
        test_dict = ui._compare(model=test_model, comparison=test_data)
        ui._print_differences(differences=test_dict)


if __name__ == "__main__":
    unittest.main(verbosity=2)
