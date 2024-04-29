#!/usr/bin/env python3
"""Unitest for build of summary

This modules checks functions responsable for the creation of txt, csv and xlsx
files.
"""

import tempfile
from pathlib import Path
from unittest import TestCase, main

import numpy
import pandas
from cobra import Model
from cobramod.core.summary import DataModel, summary
from cobramod.debug import change_to_debug
from cobramod.test import textbook

change_to_debug()


class TestSummary(TestCase):
    def test_from_model(self):
        # Case 1 example cobra model
        data_model = DataModel.from_model(textbook)
        self.assertIsInstance(obj=data_model, cls=DataModel)

        # Case 2 Empty cobra model
        model = Model("0")
        data_model = DataModel.from_model(model=model)
        self.assertIsInstance(obj=data_model, cls=DataModel)

        for variable in vars(data_model).values():
            self.assertEqual(variable, [])

    def test_diff(self):
        # Preparation
        data_model = DataModel.from_model(textbook)
        data_model_empty = DataModel(
            {
                "reactions": [],
                "metabolites": [],
                "demands": [],
                "exchanges": [],
                "genes": [],
                "groups": [],
                "sinks": [],
            }
        )
        tiny_base = DataModel(
            {
                "reactions": ["ACALD", "ACONTa"],
                "metabolites": ["13dpg_c", "acald_c"],
                "demands": ["h2o_e", "lac__D_e"],
                "exchanges": ["EX_glc__D_e", "EX_lac__D_e"],
                "genes": ["b1241", "b0474"],
                "groups": ["PWY-1187"],
                "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
            }
        )

        tiny_changed = DataModel(
            {
                "reactions": ["ACALD", "ACKr"],
                "metabolites": ["acald_c", "2pg_c"],
                "demands": ["h2o_e", "h_e"],
                "exchanges": ["EX_lac__D_e", "EX_h2o_e"],
                "genes": ["b1241", "s0001"],
                "groups": ["PWY-1187-b"],
                "sinks": ["SK_GLUTATHIONE_c", "SK_AMMONIUM_c"],
            }
        )

        tiny_changes = DataModel(
            {
                "reactions": ["ACONTa", "ACKr"],
                "metabolites": ["13dpg_c", "2pg_c"],
                "demands": ["lac__D_e", "h_e"],
                "exchanges": ["EX_glc__D_e", "EX_h2o_e"],
                "genes": ["b0474", "s0001"],
                "groups": ["PWY-1187", "PWY-1187-b"],
                "sinks": ["SK_GLY_c", "SK_AMMONIUM_c"],
            }
        )

        # CASE 1: no differences
        diff = data_model.diff(data_model)

        for variable in vars(diff).values():
            self.assertEqual(variable, [])

        # Case 2: compare with empty models
        diff_left = data_model.diff(data_model_empty)
        diff_right = data_model_empty.diff(data_model)
        dict_left = vars(diff_left)
        dict_right = vars(diff_right)

        for key, value in vars(data_model).items():
            self.assertEqual(value, dict_left[key])
            self.assertEqual(value, dict_right[key])

        # Case 3: check with tiny_base for unidirectional similarity
        diff_left = data_model.diff(tiny_base)
        diff_right = tiny_base.diff(data_model)
        dict_left = vars(diff_left)
        dict_right = vars(diff_right)

        for key, value in vars(data_model).items():
            self.assertNotEqual(value, dict_left[key])
            self.assertNotEqual(value, dict_right[key])
            self.assertEqual(dict_left[key], dict_right[key])

        # Case 4: check correct difference identification
        diff_left = tiny_base.diff(tiny_changed)
        diff_right = tiny_changed.diff(tiny_base)
        dict_left = vars(diff_left)
        dict_right = vars(diff_right)

        for key, value in vars(tiny_changes).items():
            self.assertCountEqual(value, dict_left[key])
            self.assertCountEqual(value, dict_right[key])

    def test_to_dataframe(self):
        # Preparation
        tiny_base = DataModel(
            {
                "reactions": ["ACALD", "ACONTa"],
                "metabolites": ["13dpg_c", "acald_c"],
                "demands": ["h2o_e", "lac__D_e"],
                "exchanges": ["EX_glc__D_e", "EX_lac__D_e"],
                "genes": ["b1241", "b0474"],
                "groups": ["PWY-1187"],
                "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
            }
        )

        expected_dataframe = pandas.DataFrame(
            {
                "Model identifier": ["e_coli_core", numpy.nan],
                "Model name": ["None", numpy.nan],
                "Reactions": ["ACALD", "ACONTa"],
                "Exchange": ["EX_glc__D_e", "EX_lac__D_e"],
                "Demand": ["h2o_e", "lac__D_e"],
                "Sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
                "Metabolites": ["13dpg_c", "acald_c"],
                "Genes": ["b1241", "b0474"],
                "Groups": ["PWY-1187", numpy.nan],
                "New in Reactions": ["ACALD", "ACONTa"],
                "New in Exchange": ["EX_glc__D_e", "EX_lac__D_e"],
                "New in Demand": ["h2o_e", "lac__D_e"],
                "New in Sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
                "New in Metabolites": ["13dpg_c", "acald_c"],
                "New in Genes": ["b1241", "b0474"],
                "New in Groups": ["PWY-1187", numpy.nan],
            },
            dtype=pandas.StringDtype(),
        )

        dataframe = tiny_base._to_dataframe(textbook, tiny_base)

        # Case 1: check that the dataframe contains the correct values

        pandas.testing.assert_frame_equal(expected_dataframe, dataframe)

        # Case 2: same as before but now with deletions

        expected_dataframe = pandas.DataFrame(
            {
                "Model identifier": ["e_coli_core", numpy.nan],
                "Model name": ["None", numpy.nan],
                "Reactions": ["ACALD", "ACONTa"],
                "Exchange": ["EX_glc__D_e", "EX_lac__D_e"],
                "Demand": ["h2o_e", "lac__D_e"],
                "Sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
                "Metabolites": ["13dpg_c", "acald_c"],
                "Genes": ["b1241", "b0474"],
                "Groups": ["PWY-1187", numpy.nan],
                "New in Reactions": ["ACALD", "ACONTa"],
                "New in Exchange": ["EX_glc__D_e", "EX_lac__D_e"],
                "New in Demand": ["h2o_e", "lac__D_e"],
                "New in Sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
                "New in Metabolites": ["13dpg_c", "acald_c"],
                "New in Genes": ["b1241", "b0474"],
                "New in Groups": ["PWY-1187", numpy.nan],
                "Removed in Reactions": ["ACALD", "ACONTa"],
                "Removed in Exchange": ["EX_glc__D_e", "EX_lac__D_e"],
                "Removed in Demand": ["h2o_e", "lac__D_e"],
                "Removed in Sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
                "Removed in Metabolites": ["13dpg_c", "acald_c"],
                "Removed in Genes": ["b1241", "b0474"],
                "Removed in Groups": ["PWY-1187", numpy.nan],
            },
            dtype=pandas.StringDtype(),
        )

        dataframe = tiny_base._to_dataframe(textbook, tiny_base, tiny_base)

        pandas.testing.assert_frame_equal(expected_dataframe, dataframe)

    def test_to_excl(self):
        tiny_base = DataModel(
            {
                "reactions": ["ACALD", "ACONTa"],
                "metabolites": ["13dpg_c", "acald_c"],
                "demands": ["h2o_e", "lac__D_e"],
                "exchanges": ["EX_glc__D_e", "EX_lac__D_e"],
                "genes": ["b1241", "b0474"],
                "groups": ["PWY-1187"],
                "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
            }
        )

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "summary.xlsx"
            tiny_base.to_excl(filename, textbook, tiny_base, tiny_base)

            self.assertTrue(expr=filename.exists())

    def test_to_csv(self):
        tiny_base = DataModel(
            {
                "reactions": ["ACALD", "ACONTa"],
                "metabolites": ["13dpg_c", "acald_c"],
                "demands": ["h2o_e", "lac__D_e"],
                "exchanges": ["EX_glc__D_e", "EX_lac__D_e"],
                "genes": ["b1241", "b0474"],
                "groups": ["PWY-1187"],
                "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
            }
        )

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "summary.csv"
            tiny_base.to_csv(filename, textbook, tiny_base, tiny_base)

            self.assertTrue(expr=filename.exists())

    def test_to_txt(self):
        tiny_base = DataModel(
            {
                "reactions": ["ACALD", "ACONTa"],
                "metabolites": ["13dpg_c", "acald_c"],
                "demands": ["h2o_e", "lac__D_e"],
                "exchanges": ["EX_glc__D_e", "EX_lac__D_e"],
                "genes": ["b1241", "b0474"],
                "groups": ["PWY-1187"],
                "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
            }
        )

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "summary.txt"
            tiny_base.to_txt(filename, textbook, tiny_base, tiny_base)

            self.assertTrue(expr=filename.exists())

    def test_summary(self):
        # Preparation
        test_model = textbook.copy()
        old_values = DataModel.from_model(textbook)

        # Case 1: None as filename
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            filename = None
            summary(test_model, old_values, filename=filename)

            for _ in directory.iterdir():
                self.fail("Summary created a file although None was passed")

        # Case 2: Excel as format
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            filename = directory / "summary.xlsx"
            summary(test_model, old_values, filename=filename)

            self.assertTrue(expr=filename.exists())

            for file in directory.iterdir():
                if file != filename:
                    self.fail("Summary created additional files")

        # Case 3: csv as format
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            filename = directory / "summary.csv"
            summary(test_model, old_values, filename=filename)

            self.assertTrue(expr=filename.exists())

            for file in directory.iterdir():
                if file != filename:
                    self.fail("Summary created additional files")

        # Case 4: txt as format
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            filename = directory / "summary.txt"
            summary(test_model, old_values, filename=filename)

            self.assertTrue(expr=filename.exists())

            for file in directory.iterdir():
                if file != filename:
                    self.fail("Summary created additional files")

        # Case 5: unknown format
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            filename = directory / "summary.unknown"
            summary(test_model, old_values, filename=filename)

            for _ in directory.iterdir():
                self.fail("Summary created a file!")


if __name__ == "__main__":
    main(verbosity=2)
