#!/usr/bin/env python3
from itertools import chain
from logging import DEBUG
from pathlib import Path
import unittest

from cobramod.debug import debug_log
from cobramod import get_data
import cobramod.graph as gr


debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")

if not dir_data.exists():
    dir_data.mkdir(parents=True)


class GraphTesting(unittest.TestCase):
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
        """" Start vertex are always keys """
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
                replacement_dict=test_replacement,
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
