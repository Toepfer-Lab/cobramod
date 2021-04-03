import tempfile
from filecmp import cmp
from logging import DEBUG
from pathlib import Path
from unittest import TestCase

import numpy
import pandas

from cobramod.core.summary import DataModel
from cobramod.debug import debug_log
from cobramod.test import textbook

debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")

class TestDataModel(TestCase):
    def test_from_model(self):
        data_model = DataModel.from_model(textbook)
        self.assertIsInstance(obj=data_model, cls=DataModel)

    def test_diff(self):
        # Preparation
        data_model = DataModel.from_model(textbook)
        data_model_empty = DataModel({
            "reactions": [],
            "metabolites": [],
            "demands": [],
            "exchanges": [],
            "genes": [],
            "groups": [],
            "sinks": [],
        })
        tiny_base = DataModel({
            "reactions": ['ACALD', 'ACONTa'],
            "metabolites": ['13dpg_c', 'acald_c'],
            "demands": ['h2o_e', 'lac__D_e'],
            "exchanges": ['EX_glc__D_e', 'EX_lac__D_e'],
            "genes": ['b1241', 'b0474'],
            "groups": ["PWY-1187"],
            "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
        })

        tiny_changed = DataModel({
            "reactions": ['ACALD', 'ACKr'],
            "metabolites": ['acald_c', '2pg_c'],
            "demands": ['h2o_e', 'h_e'],
            "exchanges": ['EX_lac__D_e', 'EX_h2o_e'],
            "genes": ['b1241', 's0001'],
            "groups": ["PWY-1187-b"],
            "sinks": ["SK_GLUTATHIONE_c", "SK_AMMONIUM_c"],
        })

        tiny_changes = DataModel({
            "reactions": ['ACONTa', 'ACKr'],
            "metabolites": ['13dpg_c', '2pg_c'],
            "demands": ['lac__D_e', 'h_e'],
            "exchanges": ['EX_glc__D_e', 'EX_h2o_e'],
            "genes": ['b0474', 's0001'],
            "groups": ["PWY-1187", "PWY-1187-b"],
            "sinks": ["SK_GLY_c", "SK_AMMONIUM_c"],
        })

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
        tiny_base = DataModel({
            "reactions": ['ACALD', 'ACONTa'],
            "metabolites": ['13dpg_c', 'acald_c'],
            "demands": ['h2o_e', 'lac__D_e'],
            "exchanges": ['EX_glc__D_e', 'EX_lac__D_e'],
            "genes": ['b1241', 'b0474'],
            "groups": ["PWY-1187"],
            "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
        })

        expected_dataframe = pandas.DataFrame(
            {"Model identifier": ["e_coli_core", numpy.nan],
             "Model name": ["", numpy.nan],
             "Reactions": ['ACALD', 'ACONTa'],
             "Exchange": ['EX_glc__D_e', 'EX_lac__D_e'],
             "Demand": ['h2o_e', 'lac__D_e'],
             "Sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
             "Metabolites": ['13dpg_c', 'acald_c'],
             "Genes": ['b1241', 'b0474'],
             "Groups": ["PWY-1187", numpy.nan],
             "Changed reactions": ['ACALD', 'ACONTa'],
             "Changed exchange": ['EX_glc__D_e', 'EX_lac__D_e'],
             "Changed demand": ['h2o_e', 'lac__D_e'],
             "Changed sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
             "Changed metabolites": ['13dpg_c', 'acald_c'],
             "Changed genes": ['b1241', 'b0474'],
             "Changed groups": ["PWY-1187", numpy.nan]})

        dataframe = tiny_base._to_dataframe(textbook, tiny_base)

        # Case 1: check that the dataframe contains the correct values

        pandas.testing.assert_frame_equal(expected_dataframe, dataframe)

    def test_to_excl(self):
        tiny_base = DataModel({
            "reactions": ['ACALD', 'ACONTa'],
            "metabolites": ['13dpg_c', 'acald_c'],
            "demands": ['h2o_e', 'lac__D_e'],
            "exchanges": ['EX_glc__D_e', 'EX_lac__D_e'],
            "genes": ['b1241', 'b0474'],
            "groups": ["PWY-1187"],
            "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
        })

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "summary.xlsx"
            tiny_base.to_excl(filename, textbook, tiny_base)
            self.assertTrue(expr=filename.exists())

    def test_to_csv(self):
        tiny_base = DataModel({
            "reactions": ['ACALD', 'ACONTa'],
            "metabolites": ['13dpg_c', 'acald_c'],
            "demands": ['h2o_e', 'lac__D_e'],
            "exchanges": ['EX_glc__D_e', 'EX_lac__D_e'],
            "genes": ['b1241', 'b0474'],
            "groups": ["PWY-1187"],
            "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
        })

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "summary.csv"
            tiny_base.to_csv(filename, textbook, tiny_base)

            expected_csv = pandas.read_csv(dir_input / "test_summary.csv")
            written_csv = pandas.read_csv(filename)
            pandas.testing.assert_frame_equal(expected_csv,written_csv)

    def test_to_txt(self):
        tiny_base = DataModel({
            "reactions": ['ACALD', 'ACONTa'],
            "metabolites": ['13dpg_c', 'acald_c'],
            "demands": ['h2o_e', 'lac__D_e'],
            "exchanges": ['EX_glc__D_e', 'EX_lac__D_e'],
            "genes": ['b1241', 'b0474'],
            "groups": ["PWY-1187"],
            "sinks": ["SK_GLUTATHIONE_c", "SK_GLY_c"],
        })

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "summary.txt"
            tiny_base.to_txt(filename, textbook, tiny_base)

            self.assertTrue(expr=cmp(filename, dir_input / "test_summary.txt", shallow=False))
