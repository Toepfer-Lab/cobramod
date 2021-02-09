#!/usr/bin/env python3
"""Unit test for sub-package visualization

This module includes two TestCases:

- TestItems: Creation and behaviour of JSON objects for the Escher-schema
- TestJsonDictionary: Testing the methods inside the JsonDictionary
"""
from contextlib import suppress
from pathlib import Path
from time import sleep
from unittest import TestCase, main, skip

from escher import Builder

from cobramod.error import FoundInPairError
from cobramod.visualization.pair import PairDictionary
from cobramod.visualization.items import Node, Segment, Reaction
from cobramod.visualization.converter import JsonDictionary, _convert_string
import cobramod.visualization.mapping as mp


dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")


class TestItems(TestCase):
    """
    Behaviour of JSON objects
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
        test_string = (
            "C00001_c + 2 C00002_c --> C00009_c + C00080_c + G11113_c"
        )
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
        # CASE 0: Check instance behaviour.
        test_dict = Node(node_type="midmarker", x=1, y=2)
        test_dict_2 = Node(node_type="midmarker", x=1, y=2)
        self.assertIsNot(expr1=test_dict, expr2=test_dict_2)
        # CASE 1a: Regular marker
        test_dict = Node(node_type="midmarker", x=1, y=2)
        self.assertNotIn(member="x_label", container=test_dict.keys())
        self.assertEqual(first=test_dict["node_type"], second="midmarker")
        # CASE 1b: Regular metabolite
        test_dict = Node(
            node_type="metabolite", x=1, y=2, label_y=1, label_x=2
        )
        self.assertIn(member="label_x", container=test_dict.keys())
        self.assertEqual(first=test_dict["node_type"], second="metabolite")

    def test_Segment(self):
        # CASE 0: Check instance behaviour.
        test_dict = Segment(from_node_id="1", to_node_id="2")
        test_dict_2 = Segment(from_node_id="1", to_node_id="2")
        self.assertIsNot(expr1=test_dict, expr2=test_dict_2)
        # CASE 1a: Regular creation of Segment. b1 and b2 are not specified
        test_dict = Segment(from_node_id="1", to_node_id="2")
        self.assertIn(member="from_node_id", container=test_dict.keys())
        self.assertIn(member="b2", container=test_dict.keys())
        # CASE 1b: b1 and b2 are specified
        test_dict = Segment(
            from_node_id="1",
            to_node_id="2",
            b1={"x": 10, "y": 100},
            b2={"x": 20, "y": 200},
        )
        self.assertIn(member="x", container=test_dict["b1"].keys())
        self.assertIn(member="y", container=test_dict["b2"].keys())

    def test_Reaction(self):
        # CASE 0: Check instance behaviour.
        test_dict = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
        )
        test_dict_2 = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=False,
            label_x=10,
            label_y=20,
        )
        self.assertIsNot(expr1=test_dict, expr2=test_dict_2)
        # CASE 1a: Regular creation.
        test_dict = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
        )
        self.assertIn(member="metabolites", container=test_dict.keys())
        # CASE 1b: Regular creation with a Segment
        test_dict = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
            segments={"0": Segment(from_node_id="0", to_node_id="1")},
        )
        self.assertIn(member="b1", container=test_dict["segments"]["0"].keys())
        # CASE 2: method add_metabolite
        test_dict.add_metabolite(bigg_id="test_metabolite", coefficient=1)
        self.assertEqual(
            first="test_metabolite",
            second=test_dict["metabolites"][0]["bigg_id"],
        )


class TestJsonDictionary(TestCase):
    """
    Methods for the JsonDictionary.
    """

    def test___init__(self):
        # CASE 0: Checking behaviour with two instances
        test_dict = JsonDictionary()
        test_dict_2 = JsonDictionary()
        test_dict["reactions"]["0"] = "test_string"
        self.assertIsNot(expr1=test_dict, expr2=test_dict_2)
        self.assertRaises(KeyError, lambda: test_dict_2["reactions"]["0"])
        # CASE 1: creation of dictionary without extra arguments.
        test_dict = JsonDictionary()
        self.assertEqual(first={}, second=test_dict["reactions"])
        self.assertEqual(first="", second=test_dict["head"]["map_name"])
        self.assertEqual(first=1500, second=test_dict["canvas"]["width"])
        # CASE 2: creation of dictionary passing arguments
        test_dict = JsonDictionary(
            canvas={"x": 0, "y": 0, "width": 2000, "height": 2000}
        )
        self.assertEqual(first=2000, second=test_dict["canvas"]["width"])
        self.assertEqual(first=2000, second=test_dict["canvas"]["height"])

    def test__get_last_number(self):
        # CASE 1: get last number in JsonDictionary. Reactions are not included
        # but their segments
        test_dict = JsonDictionary()
        test_dict["reactions"]["99"] = Reaction(
            name="test_reaction",
            bigg_id="test_identifier",
            reversibility=True,
            label_x=100,
            label_y=200,
            segments={"0": Segment(from_node_id="0", to_node_id="1")},
        )
        test_dict["nodes"]["1"] = Node(node_type="midmarker", x=1, y=2)
        self.assertEqual(
            first=2, second=test_dict._get_last_number(item="nodes")
        )

    def test_create_reaction(self):
        # CASE 1a: method create_reaction. Regular.
        test_dict = JsonDictionary()
        test_reaction = test_dict.create_reaction(
            name="test_reaction",
            identifier="test_identifier",
            reversibility=True,
            segments=dict(),
        )
        self.assertEqual(
            first="test_identifier", second=test_reaction["bigg_id"]
        )

    def test_add_blank(self):
        # CASE 1: Regular initialization.
        test_dict = JsonDictionary()
        test_dict.add_blank()
        self.assertEqual(first=len(test_dict["reactions"]), second=1)
        self.assertEqual(
            first=len(test_dict["reactions"]["0"]["metabolites"]), second=0
        )
        # Case 2: Stacking empty spaces.
        test_dict.add_blank()
        test_dict.add_blank()
        self.assertEqual(first=len(test_dict["reactions"]), second=3)
        self.assertEqual(
            first=len(test_dict["reactions"]["1"]["metabolites"]), second=0
        )
        self.assertEqual(
            first=len(test_dict["reactions"]["2"]["metabolites"]), second=0
        )

    def test_add_reaction(self):
        # CASE 1: test that edges are working properly
        test_dict = JsonDictionary(
            canvas={"x": 0, "y": 0, "width": 1000, "height": 600}
        )
        test_dict.add_reaction(
            string="C00002_c + C00009_c --> C00227_c + C00003_c",
            identifier="Reaction-A",
        )
        test_dict.add_reaction(
            string="c00003_c --> C00228_c", identifier="Reaction-B"
        )
        test_dict.add_reaction(
            string="C00009_c + C00228_c--> C00004_c", identifier="Reaction-C"
        )
        test_dict.add_reaction(
            string="C00004_c + C00011_c --> C00001_c + C00200_c",
            identifier="Reaction-D",
        )
        # Reaction "1" after half of canvas, "2" before half of canvas
        self.assertGreater(a=test_dict["reactions"]["1"]["label_x"], b=500)
        self.assertLess(a=test_dict["reactions"]["2"]["label_x"], b=500)
        # Node "11" starts the third reaction (labeled "2"), since there is a
        # shared metabolite
        for segment in test_dict["reactions"]["2"]["segments"].values():
            self.assertGreater(a=int(segment["from_node_id"]), b=10)
            self.assertGreater(a=int(segment["to_node_id"]), b=10)
        # CASE 2: Catch Warning if next reaction will be out of canvas
        self.assertWarns(
            UserWarning,
            test_dict.add_reaction,
            string="2 C00200_c --> 4 C00021_c",
            identifier="Reaction-E",
        )
        # CASE 3: Test multiple reactions with different participants
        test_dict = JsonDictionary()
        test_dict.add_reaction(
            string="C00002_c + C00009_c --> C00227_c + C00003_c",
            identifier="Reaction-A",
        )
        test_dict.add_reaction(
            string="C00003_c --> C00228_c", identifier="Reaction-B"
        )
        test_dict.add_reaction(
            string="C00009_c + C00228_c--> C00004_c", identifier="Reaction-C"
        )
        test_dict.add_reaction(
            string="C00004_c + C00011_c --> C00001_c + C00200_c",
            identifier="Reaction-D",
        )
        test_dict.add_reaction(
            string="2 C00200_c --> 4 C00021_c", identifier="Reaction-E"
        )
        test_dict.add_reaction(
            string="2 C00021_c + C00002_c--> C00033_c", identifier="Reaction-F"
        )
        test_dict.add_reaction(
            string="4 C00228_c + C00033_c + C00009_c --> C00011_c + "
            + "2 C00034_c + C00004_c + C00226_c",
            identifier="Reaction-G",
        )

    def test_reaction_and_blank(self):
        # CASE 1: Mixing both, one after one
        test_dict = JsonDictionary()
        test_dict.add_reaction(
            string="C00002_c + C00009_c --> C00227_c + C00003_c",
            identifier="Reaction-A",
        )
        test_dict.add_blank()
        self.assertEqual(first=len(test_dict["reactions"]), second=2)
        self.assertEqual(
            first=len(test_dict["reactions"]["1"]["metabolites"]), second=0
        )
        # CASE 2: stacking reactions
        test_dict.add_reaction(
            string="c00003_c --> C00228_c", identifier="Reaction-B"
        )
        test_dict.add_reaction(
            string="C00009_c + C00228_c--> C00004_c", identifier="Reaction-C"
        )
        test_dict.add_reaction(
            string="C00004_c + C00011_c --> C00001_c + C00200_c",
            identifier="Reaction-D",
        )
        test_dict.add_blank()
        self.assertEqual(first=len(test_dict["reactions"]), second=6)
        test_dict.add_reaction(
            string="2 C00200_c --> 4 C00021_c", identifier="Reaction-E"
        )
        test_dict.add_reaction(
            string="2 C00021_c + C00002_c--> C00033_c", identifier="Reaction-F"
        )
        test_dict.add_blank()
        test_dict.add_reaction(
            string="4 C00228_c + C00033_c + C00009_c --> C00011_c + "
            + "2 C00034_c + C00004_c + C00226_c",
            identifier="Reaction-G",
        )
        self.assertEqual(first=len(test_dict["reactions"]), second=10)

    def test_json_dump(self):
        # CASE 1: Simple HTML and JSON with 4 reactions
        test_dict = JsonDictionary()
        # Escher builder
        test_builder = Builder()
        test_dict.add_reaction(
            string="C00004_c + C00011_c --> C00001_c + C00200_c",
            identifier="Reaction-D",
        )
        test_dict.add_reaction(
            string="2 C00200_c --> 4 C00021_c", identifier="Reaction-E"
        )
        test_dict.add_reaction(
            string="2 C00021_c + C00002_c--> C00033_c", identifier="Reaction-F"
        )
        test_dict.add_reaction(
            string="4 C00228_c + C00033_c + C00009_c --> C00011_c + "
            + "2 C00034_c + C00004_c + C00226_c",
            identifier="Reaction-G",
        )
        # Writing the JSON
        test_string = test_dict.json_dump(indent=4)
        # Load the JSON and save the builder. Remove previous files.
        test_builder.map_json = test_string
        test_path = Path.cwd().joinpath("test_map.html")
        with suppress(FileNotFoundError):
            test_path.unlink()
        test_builder.save_html(str(test_path))
        self.assertTrue(expr=test_path.exists())

    @skip("")
    def test_visualize(self):
        test_path = Path.cwd().joinpath("test_map.html")
        # CASE 1: regular visualization without data
        test_dict = JsonDictionary()
        with suppress(FileNotFoundError):
            test_path.unlink()
        # Escher builder
        test_dict.add_reaction(
            string="C00004_c + C00011_c --> C00001_c + C00200_c",
            identifier="Reaction-D",
        )
        test_dict.add_reaction(
            string="2 C00200_c --> 4 C00021_c", identifier="Reaction-E"
        )
        test_builder = test_dict.visualize(filepath=test_path)
        sleep(1)
        self.assertEqual(first=test_builder.reaction_data, second=None)
        self.assertTrue(expr=test_path.exists())
        # CASE 2: visualization with Data
        test_dict = JsonDictionary()
        with suppress(FileNotFoundError):
            test_path.unlink()
        test_flux = {"Reaction-D": 2, "Reaction-E": -1}
        # Escher builder
        test_dict.add_reaction(
            string="C00004_c + C00011_c --> C00001_c + C00200_c",
            identifier="Reaction-D",
        )
        test_dict.add_reaction(
            string="2 C00200_c --> 4 C00021_c", identifier="Reaction-E"
        )
        test_dict.reaction_data = test_flux
        test_builder = test_dict.visualize(filepath=test_path)
        sleep(1)
        self.assertEqual(
            first=test_builder.reaction_data["Reaction-D"], second=2
        )
        self.assertTrue(expr=test_path.exists())
        # CASE 3: Check if blanks appears in visualization.
        test_dict = JsonDictionary()
        with suppress(FileNotFoundError):
            test_path.unlink()
        # Escher builder
        test_dict.add_reaction(
            string="C00002_c + C00009_c --> C00227_c + C00003_c",
            identifier="Reaction-A",
        )
        test_dict.add_blank()
        test_dict.add_reaction(
            string="c00003_c --> C00228_c", identifier="Reaction-B"
        )
        test_dict.add_reaction(
            string="C00009_c + C00228_c--> C00004_c", identifier="Reaction-C"
        )
        test_dict.add_reaction(
            string="C00004_c + C00011_c --> C00001_c + C00200_c",
            identifier="Reaction-D",
        )
        test_dict.add_blank()
        test_dict.add_reaction(
            string="2 C00200_c --> 4 C00021_c", identifier="Reaction-E"
        )
        test_dict.add_reaction(
            string="2 C00021_c + C00002_c--> C00033_c", identifier="Reaction-F"
        )
        test_dict.add_blank()
        test_dict.add_reaction(
            string="4 C00228_c + C00033_c + C00009_c --> C00011_c + "
            + "2 C00034_c + C00004_c + C00226_c",
            identifier="Reaction-G",
        )
        test_builder = test_dict.visualize(filepath=test_path)
        sleep(1)
        # Blanks are removed in visualization
        self.assertEqual(first=len(test_dict["reactions"]), second=7)
        # CASE 4: Blanks and information
        test_flux = {"Reaction-A": 2, "Reaction-D": 2, "Reaction-E": -1}
        test_dict.reaction_data = test_flux
        test_builder = test_dict.visualize(filepath=test_path)

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
        # CASE 1: Simple Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_matrix = mp.unformatted_matrix(graph=test_dict)
        self.assertIn(member=["R1", "R2", "R3"], container=test_matrix)
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


if __name__ == "__main__":
    main(verbosity=2)
