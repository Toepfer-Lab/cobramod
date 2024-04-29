#!/usr/bin/env python3
"""Unittest for module graph

In this modules, the new algorithm is tested. The TestCase GraphTesting checks
that the behavior of all important functions works as intended
"""

import unittest
from pathlib import Path
from typing import Any

import cobramod.core.graph as gr
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.debug import change_to_debug
from cobramod.error import GraphKeyError

# Debug must be set in level DEBUG for the test
change_to_debug()
# Setting directory for data
dir_data = Path(__file__).resolve().parent.joinpath("data")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class GraphTesting(unittest.TestCase):
    """
    This new TestCase Checks that the behavior of the new algorithm works as
    intended.
    """

    def test_find_missing(self):
        # CASE 1: Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        gr.find_missing(graph=test_dict)
        # CASE 1b: Simple Lineal, Error
        test_dict = {"R1": ("R2", "R4"), "R2": "R3", "R3": None}
        self.assertRaises(GraphKeyError, gr.find_missing, graph=test_dict)

    def test_find_cycle(self):
        # CASE: Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_answer = gr.find_cycle(graph=test_dict, key="R1", visited=[])
        self.assertListEqual(test_answer, [])

        # CASE: cyclic
        test_dict: dict[str, Any] = {"R1": "R2", "R2": "R3", "R3": "R1"}
        test_answer = gr.find_cycle(graph=test_dict, key="R1", visited=[])
        self.assertEqual(first=["R1", "R2", "R3"], second=test_answer)

        # CASE: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": None,
            "R5": None,
        }
        test_answer = gr.find_cycle(graph=test_dict, key="R0", visited=[])
        self.assertEqual(test_answer, [])

        # CASE: Complex Cyclic
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R0",
            "R5": None,
        }
        test_answer = gr.find_cycle(graph=test_dict, key="R0", visited=[])
        self.assertEqual(
            first=["R0", "R1", "R2", "R3", "R4"], second=test_answer
        )
        # CASE: Complex full Cycle
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R0",
            "R5": "R4",
        }
        # Behavior of tuples
        test_answer = gr.find_cycle(graph=test_dict, key="R1", visited=[])
        self.assertEqual(first=len(test_answer), second=4)

    def test_return_cycles(self):
        # CASE: Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_list = gr.return_cycles(graph=test_dict)
        self.assertFalse(expr=test_list)

        # CASE: Complex Cyclic
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R0",
            "R5": None,
        }
        test_list = gr.return_cycles(graph=test_dict)
        self.assertIn(
            member=["R0", "R1", "R2", "R3", "R4"], container=test_list
        )

        # CASE: Graph with cycle, from Biocyc (GLUCONEO-PWY)
        test_dict = {
            "R1": "R13",
            "R2": "R3",
            "R3": ("R2", "R1"),
            "R4": None,
            "R5": "R4",
            "R6": "R5",
            "R7": "R6",
            "R8": "R7",
            "R9": "R8",
            "R10": None,
            "R11": "R10",
            "R12": None,
            "R13": None,
        }
        test_list = gr.return_cycles(graph=test_dict)
        self.assertIn(member=["R3", "R2"], container=test_list)
        # CASE: Complex Cycle (CALVIN-PWY)
        test_dict = {
            "R1": ("R11", "R2"),
            "R2": "R10",
            "R3": "R2",
            "R4": "R7",
            "R5": "R4",
            "R6": "R5",
            "R7": ("R8", "R13"),
            "R8": "R9",
            "R9": "R1",
            "R10": "R6",
            "R11": None,
            "R12": "R3",
            "R13": "R12",
        }
        test_list = gr.return_cycles(graph=test_dict)
        self.assertIn(
            ["R1", "R2", "R10", "R6", "R5", "R4", "R7", "R13", "R12", "R3"],
            test_list,
        )

    def test_cut_cycle(self):
        # CASE: Simple cut
        test_dict = {"R1": "R2", "R2": "R3", "R3": "R1"}
        gr.cut_cycle(graph=test_dict, key="R1")
        self.assertDictEqual(
            d1=test_dict, d2={"R1": None, "R2": "R3", "R3": "R1"}
        )

    def test_back(self):
        # CASE: Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_answer = gr.back(graph=test_dict, value="R3", path=[])
        self.assertListEqual(list1=test_answer, list2=["R1", "R2"])

        # CASE: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": None,
            "R5": None,
        }
        # For R5
        test_answer = gr.back(graph=test_dict, value="R5", path=[])
        self.assertListEqual(list1=test_answer, list2=["R0", "R1"])

        # For R4
        test_answer = gr.back(graph=test_dict, value="R4", path=[])
        self.assertListEqual(list1=test_answer, list2=["R0", "R1", "R2", "R3"])

        # CASE: Rest with stop list
        test_dict = {
            "R7": ("R8", "R11"),
            "R8": "R10",
            "R9": None,
            "R10": None,
            "R11": None,
        }
        test_answer = gr.back(
            graph=test_dict,
            value="R10",
            path=[],
            stop_list=["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],
        )
        self.assertListEqual(list1=test_answer, list2=["R7", "R8"])

        # CASE: Rest with stop list. Single element
        test_answer = gr.back(
            graph=test_dict,
            value="R9",
            path=[],
            stop_list=["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],
        )
        self.assertListEqual(list1=test_answer, list2=[])

    def test_verify_paths(self):
        # CASE: Complex Cycle (after modification)
        test_dict = {"R1": ("R11", "R2"), "R2": None, "R3": "R2", "R11": None}
        test_list = [["R1", "R2"]]
        test_answer = gr.verify_paths(paths=test_list, graph=test_dict)
        self.assertCountEqual(first=test_answer, second=["R3", "R11"])

    def test_get_paths(self):
        # CASE: Unrelated
        test_dict = {"R1": None, "R2": None, "R3": None}
        test_list = gr.get_paths(graph=test_dict, stop_list=[])
        self.assertEqual(first=len(test_list), second=3)
        self.assertIn(member=["R1"], container=test_list)
        self.assertIn(member=["R2"], container=test_list)
        self.assertIn(member=["R3"], container=test_list)

        # CASE: Simple Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_list = gr.get_paths(graph=test_dict, stop_list=[])
        self.assertListEqual(
            list1=["R1", "R2", "R3"], list2=max(test_list, key=len)
        )
        # CASE: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": None,
            "R5": None,
        }
        test_list = gr.get_paths(graph=test_dict, stop_list=[])
        self.assertListEqual(
            list1=["R0", "R1", "R2", "R3", "R4"], list2=max(test_list, key=len)
        )
        # CASE: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R7"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R5",
            "R5": ("R6", "R9"),
            "R6": "R12",
            "R7": ("R8", "R11"),
            "R8": "R10",
            "R9": None,
            "R10": None,
            "R11": None,
            "R12": None,
        }
        test_list = gr.get_paths(graph=test_dict, stop_list=[])
        self.assertEqual(first=len(test_list), second=4)
        self.assertListEqual(
            list1=["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],
            list2=max(test_list, key=len),
        )
        # CASE: Using Rest
        test_dict = {
            "R7": ("R8", "R11"),
            "R8": "R10",
            "R9": None,
            "R10": None,
            "R11": None,
        }
        stop_list = (["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],)
        test_list = gr.get_paths(graph=test_dict, stop_list=stop_list)
        self.assertEqual(first=len(test_list), second=3)
        self.assertListEqual(
            list1=["R7", "R8", "R10"], list2=max(test_list, key=len)
        )
        # CASE: Complex Cycle (after modification)
        test_dict = {
            "R1": ("R11", "R2"),
            "R2": None,
            "R3": "R2",
            "R4": "R7",
            "R5": "R4",
            "R6": "R5",
            "R7": ("R8", "R13"),
            "R8": "R9",
            "R9": "R1",
            "R10": "R6",
            "R11": None,
            "R12": "R3",
            "R13": "R12",
        }
        test_list = gr.get_paths(graph=test_dict, stop_list=stop_list)
        # TODO: due to temporal fix, list is now 5 instead of 3
        self.assertEqual(first=len(test_list), second=5)

    def test_get_mapping(self):
        # CASE 0a: Single element
        test_dict = {"R1": None}
        test_list = gr.get_mapping(graph=test_dict, stop_list=[], new=[])
        self.assertListEqual(list1=test_list, list2=[["R1"]])
        # CASE 0b: Unrelated
        test_dict = {"R1": None, "R2": None, "R3": None}
        test_list = gr.get_mapping(graph=test_dict, stop_list=[], new=[])
        self.assertIn(member=["R1"], container=test_list)
        self.assertIn(member=["R2"], container=test_list)
        self.assertIn(member=["R3"], container=test_list)
        # CASE 1: Simple Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_list = gr.get_mapping(graph=test_dict, stop_list=[], new=[])
        self.assertListEqual(list1=test_list, list2=[["R1", "R2", "R3"]])
        # CASE 2: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": None,
            "R5": None,
        }
        test_list = gr.get_mapping(graph=test_dict, stop_list=[], new=[])
        self.assertListEqual(
            list1=test_list[0], list2=["R0", "R1", "R2", "R3", "R4"]
        )
        self.assertListEqual(list1=test_list[1], list2=["R5"])
        # CASE 2b: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R7"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R5",
            "R5": ("R6", "R9"),
            "R6": "R12",
            "R7": ("R8", "R11"),
            "R8": "R10",
            "R9": None,
            "R10": None,
            "R11": None,
            "R12": None,
        }
        test_list = gr.get_mapping(graph=test_dict, stop_list=[], new=[])
        # Test the two longest maps
        self.assertListEqual(
            list1=test_list[0],
            list2=["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],
        )
        self.assertListEqual(list1=test_list[1], list2=["R7", "R8", "R10"])


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
