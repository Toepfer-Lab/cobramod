#!/usr/bin/env python3
"""Unit test for sub-package visualization

This module includes two TestCases:

- TestItems: Creation and behavior of JSON objects for the Escher-schema
- TestJsonDictionary: Testing the methods inside the JsonDictionary
"""

import unittest
from contextlib import suppress
from pathlib import Path

import cobramod.visualization.mapping as mp
from cobramod.error import FoundInPairError
from cobramod.visualization.converter import (
    JsonDictionary,
    Position,
    _convert_string,
)
from cobramod.visualization.escher import EscherIntegration
from cobramod.visualization.items import Node, Reaction, Segment
from cobramod.visualization.pair import PairDictionary

dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class TestItems(unittest.TestCase):
    """
    Behavior of JSON objects
    """

    def test__convert_string(self):
        # CASE 1: Simple reaction, reversible
        test_string = "C00002_c + C00033_c <=> C00227_c + G11113_c"
        test_dict = _convert_string(string=test_string)
        self.assertEqual(first=test_dict["C00002_c"], second=-1)
        self.assertEqual(first=test_dict["C00227_c"], second=1)
        # CASE 2: Simple reaction, irreversible
        test_string = "C00002_c + C00033_c <-- C00227_c + G11113_c"
        test_dict = _convert_string(string=test_string)
        self.assertEqual(first=test_dict["C00002_c"], second=-1)
        self.assertEqual(first=test_dict["C00227_c"], second=1)
        # CASE 3: Reaction with multiple coefficients.
        test_string = "C00001_c + 2 C00002_c --> C00009_c + C00080_c + G11113_c"
        test_dict = _convert_string(string=test_string)
        self.assertEqual(first=test_dict["C00002_c"], second=-2)
        self.assertEqual(first=test_dict["C00080_c"], second=1)

    def test_PairDictionary(self):
        # CASE 1a: Raise FoundInPairError, if common key.
        test_dict_a = PairDictionary(one=1)
        test_dict_b = PairDictionary(pair=test_dict_a)
        test_dict_b.set_pair(pair=test_dict_a)
        self.assertRaises(FoundInPairError, test_dict_b.update, one=1)
        # CASE 1b: Raise FoundInPairError, if given pair has common key.
        test_dict_a = PairDictionary(one=1)
        test_dict_b = PairDictionary(one=1)
        self.assertRaises(
            FoundInPairError, test_dict_b.set_pair, pair=test_dict_a
        )
        # CASE 2: setting a pair
        test_dict_a = PairDictionary(one=1)
        test_dict_b = PairDictionary()
        test_dict_a.set_pair(pair=test_dict_b)
        self.assertIs(expr1=test_dict_a.pair, expr2=test_dict_b)

    def test_Node(self):
        # CASE 0: Check instance behavior.
        test_class = Node(node_type="midmarker", x=1, y=2)
        test_class_2 = Node(node_type="midmarker", x=1, y=2)
        self.assertIsNot(expr1=test_class, expr2=test_class_2)
        # CASE 1a: Regular marker
        test_class = Node(node_type="midmarker", x=1, y=2)
        self.assertNotIn(member="x_label", container=test_class.keys())
        self.assertEqual(first=test_class["node_type"], second="midmarker")
        # CASE 1b: Regular metabolite
        test_class = Node(
            node_type="metabolite", x=1, y=2, label_y=1, label_x=2
        )
        self.assertIn(member="label_x", container=test_class.keys())
        self.assertEqual(first=test_class["node_type"], second="metabolite")

    def test_Segment(self):
        # CASE 0: Check instance behavior.
        test_class = Segment(from_node_id="1", to_node_id="2")
        test_class_2 = Segment(from_node_id="1", to_node_id="2")
        self.assertIsNot(expr1=test_class, expr2=test_class_2)
        # CASE 1a: Regular creation of Segment. b1 and b2 are not specified
        test_class = Segment(from_node_id="1", to_node_id="2")
        self.assertIn(member="from_node_id", container=test_class.keys())
        self.assertIn(member="b2", container=test_class.keys())
        # CASE 1b: b1 and b2 are specified
        test_class = Segment(
            from_node_id="1",
            to_node_id="2",
            b1={"x": 10, "y": 100},
            b2={"x": 20, "y": 200},
        )
        self.assertIn(member="x", container=test_class["b1"].keys())
        self.assertIn(member="y", container=test_class["b2"].keys())

    def test_Reaction(self):
        # CASE 0: Check instance behavior.
        test_class = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
        )
        test_class_2 = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=False,
            label_x=10,
            label_y=20,
        )
        self.assertIsNot(expr1=test_class, expr2=test_class_2)
        # CASE 1a: Regular creation.
        test_class = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
        )
        self.assertIn(member="metabolites", container=test_class.keys())
        # CASE 1b: Regular creation with a Segment
        test_class = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
            segments={"0": Segment(from_node_id="0", to_node_id="1")},
        )
        self.assertIn(member="b1", container=test_class["segments"]["0"].keys())
        # CASE 2: method add_metabolite
        test_class.add_metabolite(bigg_id="test_metabolite", coefficient=1)
        self.assertEqual(
            first="test_metabolite",
            second=test_class["metabolites"][0]["bigg_id"],
        )


class TestJsonDictionary(unittest.TestCase):
    """
    Methods for the JsonDictionary.
    """

    def test___init__(self):
        # CASE 0: Checking behavior with two instances
        test_class = JsonDictionary()
        test_class_2 = JsonDictionary()
        test_class["reactions"]["0"] = "test_string"
        self.assertIsNot(expr1=test_class, expr2=test_class_2)
        self.assertRaises(KeyError, lambda: test_class_2["reactions"]["0"])
        # CASE 1: creation of dictionary without extra arguments.
        test_class = JsonDictionary()
        self.assertEqual(first={}, second=test_class["reactions"])
        self.assertEqual(first="", second=test_class["head"]["map_name"])

    def test__get_last_number(self):
        # CASE 1: get last number in JsonDictionary. Reactions are not included
        # but their segments
        test_class = JsonDictionary()
        test_class["reactions"]["99"] = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
            segments={"0": Segment(from_node_id="0", to_node_id="1")},
        )
        test_class["nodes"]["1"] = Node(node_type="midmarker", x=1, y=2)
        self.assertEqual(
            first=2, second=test_class._get_last_number(item="nodes")
        )

    def test_get_column_reactions(self):
        # Preparing tests
        test_class = JsonDictionary()
        test_class._overview["R1"] = {"position": Position(row=1, column=1)}
        test_class._overview["R2"] = {"position": Position(row=2, column=1)}
        test_class._overview["R3"] = {"position": Position(row=3, column=0)}
        test_class._overview["R4"] = {"position": Position(row=4, column=1)}
        test_class._overview["R5"] = {"position": Position(row=5, column=0)}
        test_class._overview["R6"] = {"position": Position(row=2, column=0)}
        test_class._overview["R7"] = {"position": Position(row=2, column=3)}
        # CASE 1: Retrieving columns
        test_list = test_class._get_matrix_reactions(vertical=False, position=1)
        self.assertCountEqual(first=["R1", "R2", "R4"], second=test_list)
        # CASE 2: Retrieving Row
        test_list = test_class._get_matrix_reactions(vertical=True, position=2)
        self.assertCountEqual(first=["R2", "R6", "R7"], second=test_list)

    def test__get_products(self):
        # Preparing tests
        test_class = JsonDictionary()
        test_class._overview["R1"] = {"index": "0"}
        test_class._overview["R2"] = {"index": "1"}
        test_class._overview["R3"] = {"index": "2"}
        test_class._overview["R4"] = {"index": "3"}
        test_class["reactions"] = {
            "0": {
                "bigg_id": "R1",
                "metabolites": [
                    {"bigg_id": "C01001_c", "coefficient": -1},
                    {"bigg_id": "C01002_c", "coefficient": 1},
                ],
            },
            "1": {
                "bigg_id": "R2",
                "metabolites": [
                    {"bigg_id": "C01002_c", "coefficient": -1},
                    {"bigg_id": "C02001_c", "coefficient": 1},
                ],
            },
            "2": {
                "bigg_id": "R3",
                "metabolites": [
                    {"bigg_id": "C02001_c", "coefficient": -1},
                    {"bigg_id": "C03001_c", "coefficient": 1},
                    {"bigg_id": "C03002_c", "coefficient": 1},
                ],
            },
            "3": {
                "bigg_id": "R4",
                "metabolites": [
                    {"bigg_id": "C03001_c", "coefficient": -1},
                    {"bigg_id": "C04001_c", "coefficient": 1},
                    {"bigg_id": "C04002_c", "coefficient": 1},
                ],
            },
        }
        # CASE 1: Empty list
        test_dict = test_class._get_products(reactions=[])
        self.assertFalse(test_dict)
        # CASE 1: Simple reactions
        test_dict = test_class._get_products(reactions=["R1"])
        self.assertDictEqual(d1={"R1": ["C01002_c"]}, d2=test_dict)
        # CASE 1: Mix with complex reactions
        test_dict = test_class._get_products(reactions=["R1", "R2", "R4"])
        self.assertDictEqual(
            d1={
                "R1": ["C01002_c"],
                "R2": ["C02001_c"],
                "R4": ["C04001_c", "C04002_c"],
            },
            d2=test_dict,
        )

    def test__find_shared(self):
        # Preparing tests
        test_class = JsonDictionary()
        test_class._overview["R1"] = {
            "nodes": {"C01001_c": "10", "C01002_c": "11"}
        }
        test_class._overview["R2"] = {
            "nodes": {"C01002_c": "11", "C02001_c": "21"}
        }
        test_class._overview["R3"] = {
            "nodes": {"C02001_c": "21", "C03001_c": "31"}
        }
        test_dict = {"R1": ["C01002_c"], "R2": ["C02001_c"]}
        # CASE 1: Nothing is found
        # node_number, old_reaction
        test_tuple = test_class._find_shared(
            metabolite="C01003_c", products=test_dict
        )
        self.assertFalse(test_tuple[0])
        self.assertFalse(test_tuple[1])
        # CASE 2: regular finding
        # node_number, old_reaction
        test_tuple = test_class._find_shared(
            metabolite="C01002_c", products=test_dict
        )
        self.assertEqual(first=test_tuple[0], second="11")
        self.assertEqual(first=test_tuple[1], second="R1")

    def test_add_reaction(self):
        # CASE 1: Matrix 2x2, Unrelated
        test_class = JsonDictionary()
        test_class.add_reaction(
            string="C01001_c + C01002_c --> C01003_c + C01004_c",
            identifier="R1",
            name="Reaction-1",
            row=0,
            column=0,
            vertical=False,
        )
        test_class.add_reaction(
            string="C02001_c --> C02002_c",
            identifier="R2",
            name="Reaction-2",
            row=0,
            column=1,
            vertical=False,
        )
        test_class.add_reaction(
            string="C03001_c + C03002_c_c--> C03003_c",
            identifier="R3",
            name="Reaction-3",
            row=1,
            column=0,
            vertical=False,
        )
        test_class.add_reaction(
            string="C04001_c + C04002_c --> C04003_c + C04004_c",
            identifier="R4",
            name="Reaction-4",
            row=1,
            column=1,
            vertical=False,
        )
        # Reaction "R2" before two reactions boxes in x (650 * 2)
        # Reaction "R3" before one reaction box in x (650)
        self.assertLess(a=test_class["reactions"]["1"]["label_x"], b=650 * 2)
        self.assertLess(a=test_class["reactions"]["2"]["label_x"], b=650)
        # CASE 3: Different positions. Matrix 1x4
        test_class = JsonDictionary()
        test_class.add_reaction(
            string="C01001_c + C01002_c --> C01003_c + C01004_c",
            identifier="R1",
            name="Reaction-1",
            row=0,
            column=0,
            vertical=False,
        )
        test_class.add_reaction(
            string="C01003_c --> C02001_c",
            identifier="R2",
            name="Reaction-2",
            row=0,
            column=1,
            vertical=False,
        )
        test_class.add_reaction(
            string="C02001_c + C03001_c_c--> C03002_c",
            identifier="R3",
            name="Reaction-3",
            row=0,
            column=2,
            vertical=False,
        )
        test_class.add_reaction(
            string="C03002_c + C04001_c --> C04002_c + C04003_c",
            identifier="R4",
            name="Reaction-4",
            row=0,
            column=3,
            vertical=False,
        )
        # Reaction "R2" before two reactions boxes in x (550*2)
        # Reaction "R3" before three reaction boxes in x ()
        self.assertLess(a=test_class["reactions"]["1"]["label_x"], b=550 * 2)
        self.assertLess(a=test_class["reactions"]["2"]["label_x"], b=551 * 3)
        # Shared metabolite is node "2" in R1 and R2
        for segment in test_class["reactions"]["1"]["segments"].values():
            self.assertGreaterEqual(a=int(segment["from_node_id"]), b=2)
            self.assertGreaterEqual(a=int(segment["to_node_id"]), b=2)
            self.assertLessEqual(a=int(segment["to_node_id"]), b=11)
            self.assertLessEqual(a=int(segment["to_node_id"]), b=11)
        # TODO: CASE 2 connecting rows

    def test_color_grading(self):
        test_class = JsonDictionary()
        test_class.flux_solution = {"A": -2, "B": -1, "C": 0, "D": 1, "E": 2}

        expected_scale = [
            {"type": "value", "value": 1.0, "color": "rgb(237,192,110)"},
            {"type": "value", "value": 2.0, "color": "rgb(255,165,0)"},
            {"type": "value", "value": 0, "color": "rgb(220,220,220)"},
            {"type": "value", "value": -1.0, "color": "rgb(110,174,110)"},
            {"type": "value", "value": -2.0, "color": "rgb(0,128,0)"},
        ]

        # Case 1: 2 colors with colors defined as str or via rgb values
        color = ["orange", "green"]
        test_class.color_grading(color=color)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        color = [[255, 165, 0], [0, 128, 0]]
        test_class.color_grading(color=color)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        # Case 2: 2 colors with one defined as rgb and one defined as str

        color = ["orange", [0, 128, 0]]
        test_class.color_grading(color=color)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        color = [[255, 165, 0], "green"]
        test_class.color_grading(color=color)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        # Case 3: only one color defined
        # redefine expected scale
        expected_scale = [
            {"type": "value", "value": 1.0, "color": "rgb(237,192,110)"},
            {"type": "value", "value": 2.0, "color": "rgb(255,165,0)"},
            {"type": "value", "value": 0, "color": "rgb(220,220,220)"},
            {"type": "value", "value": -1.0, "color": "rgb(220,220,220)"},
            {"type": "value", "value": -2.0, "color": "rgb(220,220,220)"},
        ]

        color = ["orange", None]
        test_class.color_grading(color=color)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        # redefine expected scale
        expected_scale = [
            {"type": "value", "value": 1.0, "color": "rgb(220,220,220)"},
            {"type": "value", "value": 2.0, "color": "rgb(220,220,220)"},
            {"type": "value", "value": 0, "color": "rgb(220,220,220)"},
            {"type": "value", "value": -1.0, "color": "rgb(110,174,110)"},
            {"type": "value", "value": -2.0, "color": "rgb(0,128,0)"},
        ]

        color = ["sdfgdfvcbv", [0, 128, 0]]
        test_class.color_grading(color=color)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        # Case 4: 2 colors and min_max defined
        # redefine expected scale
        expected_scale = [
            {"type": "value", "value": 1.0, "color": "rgb(255,165,0)"},
            {"type": "value", "value": 0, "color": "rgb(220,220,220)"},
            {"type": "value", "value": -1.0, "color": "rgb(0,128,0)"},
        ]
        color = ["orange", "green"]
        test_class.color_grading(color=color, min_max=[-1, 1])

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        # Case 5: no values left due to min_max

        expected_scale = [
            {"type": "value", "value": 0, "color": "rgb(220,220,220)"}
        ]
        color = ["orange", "green"]
        test_class.color_grading(color=color, min_max=[0, 0])

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        test_class.color_grading(color=color, min_max=[5, 10])
        expected_scale = [
            {"type": "value", "value": 5, "color": "rgb(237,192,110)"},
            {"type": "value", "value": 10, "color": "rgb(255,165,0)"},
        ]

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        expected_scale = [
            {"color": "rgb(255,165,0)", "type": "value", "value": 0.1},
            {"color": "rgb(220,220,220)", "type": "value", "value": 0},
            {"color": "rgb(0,128,0)", "type": "value", "value": -0.1},
        ]
        test_class.color_grading(color=color, min_max=[-0.1, 0.1])

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        # Case 6: values close together with and without quantile
        test_class.flux_solution = {
            "A": -2,
            "B": 0.1,
            "C": 0.2,
            "D": 0.3,
            "E": 2,
        }
        expected_scale = [
            {"color": "rgb(228,206,165)", "type": "value", "value": 0.5},
            {"color": "rgb(237,192,110)", "type": "value", "value": 1.0},
            {"color": "rgb(246,178,55)", "type": "value", "value": 1.5},
            {"color": "rgb(255,165,0)", "type": "value", "value": 2.0},
            {"color": "rgb(220,220,220)", "type": "value", "value": 0},
            {"color": "rgb(0,128,0)", "type": "value", "value": -2.0},
        ]

        test_class.color_grading(color=color, quantile=False)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        expected_scale = [
            {"color": "rgb(228,206,165)", "type": "value", "value": 0.1},
            {"color": "rgb(237,192,110)", "type": "value", "value": 0.2},
            {"color": "rgb(246,178,55)", "type": "value", "value": 0.3},
            {"color": "rgb(255,165,0)", "type": "value", "value": 2.0},
            {"color": "rgb(220,220,220)", "type": "value", "value": 0},
            {"color": "rgb(0,128,0)", "type": "value", "value": -2.0},
        ]
        test_class.color_grading(color=color, quantile=True)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        # Case 7: no flux
        test_class.flux_solution = {}

        expected_scale = [
            {"color": "rgb(220,220,220)", "type": "value", "value": 0}
        ]

        test_class.color_grading(color=color)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

        test_class.color_grading(color=color, quantile=True)

        self.assertCountEqual(
            test_class.reaction_scale,
            expected_scale,
            msg="Calculated color scale differs from the expected one",
        )

    def test_json_dump(self):
        # CASE 1: Simple HTML and JSON with 4 reactions
        test_class = JsonDictionary()
        # Escher builder
        test_builder = EscherIntegration()
        # Matrix 1x4
        test_class.add_reaction(
            string="C01001_c + C01002_c --> C01003_c + C01004_c",
            identifier="R1",
            name="Reaction-1",
            row=0,
            column=0,
            vertical=False,
        )
        test_class.add_reaction(
            string="C01003_c --> C02001_c",
            identifier="R2",
            name="Reaction-2",
            row=0,
            column=1,
            vertical=False,
        )
        test_class.add_reaction(
            string="C02001_c + C03001_c_c--> C03002_c",
            identifier="R3",
            name="Reaction-3",
            row=0,
            column=2,
            vertical=False,
        )
        test_class.add_reaction(
            string="C03002_c + C04001_c --> C04002_c + C04003_c",
            identifier="R4",
            name="Reaction-4",
            row=0,
            column=3,
            vertical=False,
        )
        # Writing the JSON
        test_string = test_class.json_dump(indent=4)
        # Load the JSON and save the builder. Remove previous files.
        test_builder.map_json = test_string
        test_path = Path.cwd().joinpath("test_map.html")
        with suppress(FileNotFoundError):
            test_path.unlink()
        test_builder.save_html(str(test_path))
        # ToDo Check that html actually displays something
        self.assertTrue(expr=test_path.exists())

    def test_visualize(self):
        # NOTE: visual tests
        # Settings
        test_path = Path.cwd().joinpath("test_map.html")
        test_reactions = {
            "R1": "C01001_c + C01002_c --> C01003_c + C01004_c",
            "R2": "C01003_c --> C02001_c",
            "R3": "C02001_c + C03001_c--> C03002_c",
            "R4": "C03002_c + C04001_c --> C04002_c + C04003_c",
        }
        # CASE 1: Simple single reaction
        test_class = JsonDictionary()
        test_class.graph = {"R1": None}
        test_class.reaction_strings = test_reactions
        with suppress(FileNotFoundError):
            test_path.unlink()
        test_builder = test_class.visualize(
            filepath=test_path, custom_integration=True
        )
        self.assertEqual(first=test_builder.reaction_data, second=None)
        self.assertTrue(expr=test_path.exists())
        # CASE 2: visualization with Data
        test_class = JsonDictionary()
        test_class.graph = {"R1": None}
        test_class.reaction_strings = test_reactions
        with suppress(FileNotFoundError):
            test_path.unlink()
        # Setting flux
        test_flux = {"R1": 2, "R2": -1}
        test_class.flux_solution = test_flux
        test_builder = test_class.visualize(
            filepath=test_path, custom_integration=True
        )
        self.assertEqual(first=test_builder.reaction_data["R1"], second=2)
        self.assertTrue(expr=test_path.exists())


class TestMapping(unittest.TestCase):
    """
    This TestCase checks if the mapping for the visualization has a normal
    behavior
    """

    def test_child_map(self):
        test_dict = {
            "R1": "R2",
            "R2": ("R3", "R5", "R4"),
            "R3": ("R6", "R8"),
            "R4": None,
            "R5": None,
            "R6": "R7",
            "R7": "R10",
            "R8": ("R9", "R11"),
            "R9": None,
            "R10": "R14",
            "R11": ("R12", "R13"),
            "R12": None,
            "R13": None,
            "R14": None,
        }
        test_list = [
            ["R1", "R2", "R3", "R6", "R7", "R10", "R14"],
            ["R8", "R11", "R12"],
            ["R4"],
            ["R5"],
            ["R9"],
            ["R13"],
        ]
        test_answer = mp.child_map(mapping=test_list, dictionary=test_dict)
        self.assertCountEqual(first=test_answer["1"], second=["4", "5"])
        self.assertCountEqual(first=test_answer["0"], second=["2", "3", "1"])

    def test_unformatted_matrix(self):
        # CASE 0a: Unrelated items
        test_dict = {"R1": None, "R2": None, "R3": None}
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertIn(member=["R1"], container=test_matrix)
        self.assertIn(member=["R2"], container=test_matrix)
        self.assertIn(member=["R3"], container=test_matrix)
        # CASE 0b: Single unrelated item
        test_dict = {"R1": "R2", "R2": None, "R3": None}
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertIn(member=["R1", "R2"], container=test_matrix)
        self.assertIn(member=["R3"], container=test_matrix)
        # CASE 1a: Simple Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertIn(member=["R1", "R2", "R3"], container=test_matrix)
        # CASE 1b: Simple Lineal
        test_dict = {"R1": ("R2", "R4"), "R2": "R3", "R3": None, "R4": None}
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertIn(member=["R1", "R2", "R3"], container=test_matrix)
        self.assertIn(member=[0, "R4"], container=test_matrix)
        # CASE 2: Simple Cyclic
        test_dict = {"R1": "R2", "R2": "R3", "R3": "R1"}
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertCountEqual(first=["R1", "R2", "R3"], second=test_matrix[0])
        # CASE 3a: Complex Lineal
        test_dict = {
            "R1": "R2",
            "R2": ("R3", "R5", "R4"),
            "R3": ("R6", "R8"),
            "R4": None,
            "R5": None,
            "R6": "R7",
            "R7": "R10",
            "R8": ("R9", "R11"),
            "R9": None,
            "R10": "R14",
            "R11": ("R12", "R13"),
            "R12": None,
            "R13": None,
            "R14": None,
        }
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertIn(member=[0, 0, "R5"], container=test_matrix)
        self.assertIn(member=[0, 0, "R4"], container=test_matrix)
        self.assertIn(
            member=[0, 0, 0, "R8", "R11", "R12"], container=test_matrix
        )
        # CASE 3b: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R7"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R5",
            "R5": ("R6", "R9"),
            "R6": "R12",
            "R7": "R8",
            "R8": ("R10", "R11"),
            "R9": None,
            "R10": None,
            "R11": None,
            "R12": None,
        }
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertIn(member=[0, 0, 0, 0, 0, 0, "R9"], container=test_matrix)
        self.assertIn(member=[0, 0, 0, 0, "R11"], container=test_matrix)

    def test_fill_matrix(self):
        # CASE 1: Normal use
        test_matrix = [["R1", "R2", "R3", "R4", "R5"], [0, 0, "R6"], [0, "R7"]]
        mp.fill_matrix(test_matrix, length=5)
        self.assertListEqual(list1=test_matrix[1], list2=[0, 0, "R6", 0, 0])
        self.assertListEqual(list1=test_matrix[2], list2=[0, "R7", 0, 0, 0])

    def test_format_matrix(self):
        # CASE 0a: only filling
        test_matrix = [["R1", "R2", "R3", "R4", "R5"], [0, 0, "R6"]]
        test_answer = mp.format_matrix(matrix=test_matrix.copy(), max_length=5)
        test_values = set()
        for row in test_answer:
            for i in row:
                if i != 0:
                    test_values.add(i)
        self.assertCountEqual(
            first=["R1", "R2", "R3", "R4", "R5", "R6"], second=test_values
        )
        self.assertEqual(first=len(test_answer), second=2)
        self.assertListEqual(list1=test_answer[1], list2=[0, 0, "R6", 0, 0])
        # CASE 0b: Unrelated
        test_matrix = [["R1"], ["R2"], ["R3"]]
        test_answer = mp.format_matrix(matrix=test_matrix.copy(), max_length=1)
        self.assertEqual(first=test_matrix, second=test_answer)
        # CASE 1: Simple modification
        test_matrix = [
            ["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],
            [0, 0, "R7", "R8", "R10"],
            [0, 0, 0, 0, 0, 0, "R9"],
            [0, 0, 0, 0, "R11"],
        ]
        test_answer = mp.format_matrix(matrix=test_matrix.copy(), max_length=8)
        test_values = set()
        for row in test_answer:
            for i in row:
                if i != 0:
                    test_values.add(i)
        self.assertCountEqual(
            first=[
                "R0",
                "R1",
                "R2",
                "R3",
                "R4",
                "R5",
                "R6",
                "R7",
                "R8",
                "R9",
                "R10",
                "R11",
                "R12",
            ],
            second=test_values,
        )
        self.assertEqual(first=len(test_answer), second=3)
        self.assertListEqual(
            list1=test_answer[1], list2=[0, 0, "R7", "R8", "R10", 0, "R9", 0]
        )
        # CASE 2: Complex modification
        test_matrix = [
            ["R1", "R2", "R3", "R6", "R7", "R10", "R14"],
            [0, 0, 0, "R8", "R11", "R12"],
            [0, 0, "R4"],
            [0, 0, "R5"],
            [0, 0, 0, 0, "R9"],
            [0, 0, 0, 0, 0, "R13"],
        ]
        test_answer = mp.format_matrix(matrix=test_matrix.copy(), max_length=7)
        test_values = set()
        for row in test_answer:
            for i in row:
                if i != 0:
                    test_values.add(i)
        self.assertCountEqual(
            first=[
                "R1",
                "R2",
                "R3",
                "R4",
                "R5",
                "R6",
                "R7",
                "R8",
                "R9",
                "R10",
                "R11",
                "R12",
                "R13",
                "R14",
            ],
            second=test_values,
        )
        self.assertIn(
            member=[0, 0, 0, 0, "R9", "R13", 0], container=test_answer
        )

    def test_get_all_values(self):
        # CASE 1: Simple dictionary
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_set = mp.get_all_values(dictionary=test_dict, keys=["R1"])
        self.assertCountEqual(first=test_set, second={"R2"})
        # CASE 2: Simple dictionary, multiple keys
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R7"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R5",
            "R5": None,
        }
        test_set = mp.get_all_values(
            dictionary=test_dict, keys=["R1", "R5", "R3"]
        )
        self.assertCountEqual(first=test_set, second={"R2", "R7", "R4"})

    def test_transpose(self):
        # CASE 0: Simple Matrix m*m
        test_matrix = [["R1", "R2"], ["R3", "R4"]]
        test_answer = mp.transpose(matrix=test_matrix)
        self.assertEqual(first=[["R1", "R3"], ["R2", "R4"]], second=test_answer)
        # CASE 1: Different dimensions m*n
        test_matrix = [["R1", "R2", "R3", "R4"], [0, 0, "R5", 0]]
        test_answer = mp.transpose(matrix=test_matrix)
        self.assertEqual(
            first=[["R1", 0], ["R2", 0], ["R3", "R5"], ["R4", 0]],
            second=test_answer,
        )
        # CASE 2: Matrix with multiple reactions
        test_matrix = [
            ["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],
            [0, 0, "R7", "R8", "R10", 0, "R9", 0],
            [0, 0, 0, 0, "R11", 0, 0, 0],
        ]
        test_answer = mp.transpose(matrix=test_matrix)
        self.assertEqual(
            first=[
                ["R0", 0, 0],
                ["R1", 0, 0],
                ["R2", "R7", 0],
                ["R3", "R8", 0],
                ["R4", "R10", "R11"],
                ["R5", 0, 0],
                ["R6", "R9", 0],
                ["R12", 0, 0],
            ],
            second=test_answer,
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
