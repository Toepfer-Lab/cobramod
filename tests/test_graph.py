#!/usr/bin/env python3
from itertools import chain
from logging import DEBUG
from pathlib import Path
from unittest import TestCase, main

from cobramod.debug import debug_log
from cobramod.core.retrieval import get_data
import cobramod.core.graph as gr

# Debug must be set in level DEBUG for the test
debug_log.setLevel(DEBUG)
# Setting directory for data
dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class GraphTesting(TestCase):
    def test__find_end_vertex(self):
        # CASE 1: regular lineal / 2 child edges
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E"}
        self.assertIn("E", gr._find_end_vertex(vertex_dict=test_dict))
        # CASE 2: lineal with a cycle in the middle
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G",
        }
        self.assertIn("G", gr._find_end_vertex(vertex_dict=test_dict))
        # CASE 3: cyclical (regular), other cases will be examined further
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A",
        }
        self.assertRaises(
            StopIteration, next, gr._find_end_vertex(vertex_dict=test_dict)
        )

    def test__find_start_vertex(self):
        # INFO:  Start vertex are always keys """
        # CASE 1: regular lineal
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E"}
        self.assertIn("A", gr._find_start_vertex(vertex_dict=test_dict))
        # CASE 2: lineal with 2 paths
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G",
        }
        self.assertIn("A", gr._find_start_vertex(vertex_dict=test_dict))
        self.assertIn("D", gr._find_start_vertex(vertex_dict=test_dict))
        # CASE 3: cyclical (regular), other cases will be examined further
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A",
        }
        self.assertRaises(
            StopIteration, next, gr._find_start_vertex(vertex_dict=test_dict)
        )

    def test__verify_return_end_vertex(self):
        # CASE 1: regular lineal
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E"}
        self.assertIn("E", gr._verify_return_end_vertex(vertex_dict=test_dict))
        # CASE 2: lineal with a cycle in the middle
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G",
        }
        self.assertIn("G", gr._verify_return_end_vertex(vertex_dict=test_dict))
        # CASE 3: cyclical (regular), other cases will be examined further
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A",
        }
        self.assertEqual(
            1, len(list(gr._verify_return_end_vertex(vertex_dict=test_dict)))
        )
        # CASE 4a: cyclical with reactions in between
        test_dict = {
            # cycle
            "A": "B",
            "B": "C",  # could also be G
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F",
        }
        self.assertEqual(
            1, len(list(gr._verify_return_end_vertex(vertex_dict=test_dict)))
        )
        # CASE 4b: cyclical with reactions in between (G instead of C)
        test_dict = {
            # cycle
            "A": "B",
            "B": "G",  # could also be C
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F",
        }
        self.assertEqual(
            1, len(list(gr._verify_return_end_vertex(vertex_dict=test_dict)))
        )
        # CASE 5: cyclycal with one extra out
        test_dict = {"A": "B", "B": "C", "C": "E", "D": "B", "E": "D"}
        self.assertEqual(
            1, len(list(gr._verify_return_end_vertex(vertex_dict=test_dict)))
        )

    def test_get_graph(self):
        # CASE 1: regular lineal
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E"}
        test_graph = gr.get_graph(vertex_dict=test_dict)
        self.assertIn(
            # Real path direction ->
            ["E", "D", "C", "B", "A"],
            [list(sequence) for sequence in test_graph],
        )
        # CASE 2: lineal with a cycle in the middle
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G",
        }
        test_graph = gr.get_graph(vertex_dict=test_dict)
        # 2 sequences
        self.assertEqual(2, len(list(test_graph)))
        # CASE 3: cyclical (reg,lar)
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A",
        }
        test_graph = gr.get_graph(vertex_dict=test_dict)
        self.assertEqual(
            6, len([x for x in set(chain.from_iterable(test_graph))])
        )
        # CASE 4a: cyclical with reactions in between
        # This can only be fix with the next step, checking that all edges
        # are mentioned at least on in the graph

        # CASE 4b: as 4a but replacing value of B
        test_dict = {
            # cycle
            "A": "B",
            "B": "G",
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F",
        }
        test_graph = gr.get_graph(vertex_dict=test_dict)
        self.assertEqual(
            8, len([x for x in set(chain.from_iterable(test_graph))])
        )
        # CASE 5: cyclycal with one extra out
        test_dict = {"A": "B", "B": "C", "C": "E", "D": "B", "E": "D"}
        test_graph = gr.get_graph(vertex_dict=test_dict)
        self.assertEqual(
            5, len([x for x in set(chain.from_iterable(test_graph))])
        )

    def test__return_verified_graph(self):
        # CASE 0: regular lineal graph
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E"}
        test_graph = list(
            gr._return_verified_graph(
                vertex_dict=test_dict,
                vertex_set=set(chain.from_iterable(test_dict.items())),
            )
        )
        self.assertIn(
            ["E", "D", "C", "B", "A"],
            [list(sequence) for sequence in test_graph],
        )
        # CASE 1: cyclical (regular)
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A",
        }
        test_graph = list(
            gr._return_verified_graph(
                vertex_dict=test_dict,
                vertex_set=set(chain.from_iterable(test_dict.items())),
            )
        )
        self.assertEqual(
            6, len([x for x in set(chain.from_iterable(test_graph))])
        )
        # CASE 3: cyclical with reactions in between
        test_dict = {
            # cycle
            "A": "B",
            "B": "C",  # could also be G
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F",
        }
        test_graph = list(
            gr._return_verified_graph(
                vertex_dict=test_dict,
                vertex_set=set(chain.from_iterable(test_dict.items())),
            )
        )
        # CASE 4: using replacement dict. As Key and value
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E"}
        test_replacement = {"A": "Z", "E": "F"}
        test_graph = list(
            gr._return_verified_graph(
                vertex_dict=test_dict,
                vertex_set=set(chain.from_iterable(test_dict.items())),
                replacement=test_replacement,
            )
        )
        self.assertIn("Z", set(chain.from_iterable(test_graph)))
        self.assertIn("F", set(chain.from_iterable(test_graph)))

    def test_return_from_dict(self):
        # CASE 1: regular biocyc
        test_dict = get_data(
            identifier="PWY-1187", directory=dir_data, database="META"
        )
        test_list = gr.return_graph_from_dict(data_dict=test_dict)
        # 3 end-metabolites + 1 in-between path (cyclic)
        self.assertEqual(
            first=len(set(chain.from_iterable(test_list))), second=14
        )
        # CASE 2: cycle (no repetitions)
        test_dict = get_data(
            identifier="PWY-5690", directory=dir_data, database="META"
        )
        test_list = gr.return_graph_from_dict(data_dict=test_dict)
        self.assertEqual(
            first=len(set(chain.from_iterable(test_list))), second=9
        )
        # CASE 3: Cyclical with edges in-between
        test_dict = get_data(
            identifier="CALVIN-PWY", directory=dir_data, database="META"
        )
        test_list = gr.return_graph_from_dict(data_dict=test_dict)
        # Total of 13, without repetitions
        self.assertEqual(
            first=len(set(chain.from_iterable(test_list))), second=13
        )
        # CASE 4a: regular lineal KEGG
        test_dict = get_data(
            directory=dir_data, database="KEGG", identifier="M00118"
        )
        test_list = gr.return_graph_from_dict(data_dict=test_dict)
        self.assertEqual(
            first=len(set(chain.from_iterable(test_list))), second=2
        )
        # CASE 4b: lineal with extra reactions, KEGG
        test_dict = get_data(
            directory=dir_data, database="KEGG", identifier="M00001"
        )
        test_list = gr.return_graph_from_dict(data_dict=test_dict)
        self.assertEqual(
            first=len(set(chain.from_iterable(test_list))), second=15
        )
        self.assertGreater(a=len(test_list[0]), b=4)


class NewAlgorithm(TestCase):
    """
    This new TestCase Checks that the behaviour of the new algorithm works as
    intended.
    """

    def test_find_cycle(self):
        # CASE 1: Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_answer = gr.find_cycle(graph=test_dict, key="R1", visited=[])
        self.assertIs(expr1=False, expr2=test_answer)
        # CASE 2: cyclic
        test_dict = {"R1": "R2", "R2": "R3", "R3": "R1"}
        test_answer = gr.find_cycle(graph=test_dict, key="R1", visited=[])
        self.assertCountEqual(first=["R1", "R2", "R3"], second=test_answer)
        # CASE 2a: Complex Lineal / Missing Key
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R5": None,
        }
        self.assertRaises(
            KeyError, gr.find_cycle, graph=test_dict, key="R0", visited=[]
        )
        # CASE 2a: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": None,
            "R5": None,
        }
        test_answer = gr.find_cycle(graph=test_dict, key="R0", visited=[])
        self.assertIs(expr1=False, expr2=test_answer)
        # CASE 3: Complex Cyclic
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": "R0",
            "R5": None,
        }
        test_answer = gr.find_cycle(graph=test_dict, key="R0", visited=[])
        self.assertCountEqual(
            first=["R0", "R1", "R2", "R3", "R4"], second=test_answer
        )

    def test_cut_cycle(self):
        # CASE 1: Simple cut
        test_dict = {"R1": "R2", "R2": "R3", "R3": "R1"}
        gr.cut_cycle(graph=test_dict, key="R1")
        self.assertDictEqual(
            d1=test_dict, d2={"R1": None, "R2": "R3", "R3": "R1"}
        )

    def test_back(self):
        # CASE 1: Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_answer = gr.back(graph=test_dict, value="R3", path=[])
        self.assertListEqual(list1=test_answer, list2=["R1", "R2"])
        # CASE 2: Complex Lineal
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

    def test_items(self):
        # CASE 1: Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_set = gr.items(graph=test_dict)
        self.assertCountEqual(first={"R1", "R2", "R3"}, second=test_set)
        # CASE 2: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": None,
            "R5": None,
        }
        test_set = gr.items(graph=test_dict)
        self.assertCountEqual(
            first={"R0", "R1", "R2", "R3", "R4", "R5"}, second=test_set
        )

    def test_longest_path(self):
        # CASE 1: Simple Lineal
        test_dict = {"R1": "R2", "R2": "R3", "R3": None}
        test_list = gr.longest_path(graph=test_dict)
        self.assertListEqual(list1=["R1", "R2", "R3"], list2=test_list)
        # CASE 2a: Complex Lineal
        test_dict = {
            "R0": "R1",
            "R1": ("R2", "R5"),
            "R2": "R3",
            "R3": "R4",
            "R4": None,
            "R5": None,
        }
        test_list = gr.longest_path(graph=test_dict)
        self.assertListEqual(
            list1=["R0", "R1", "R2", "R3", "R4"], list2=test_list
        )
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
        test_list = gr.longest_path(graph=test_dict)
        self.assertListEqual(
            list1=["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R12"],
            list2=test_list,
        )


if __name__ == "__main__":
    main(verbosity=2)
