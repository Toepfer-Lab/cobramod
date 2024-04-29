import json
import unittest

import cobra.core

from cobramod.test import textbook_biocyc
from cobramod.visualization.force_graph import (
    Nodes,
    GraphData,
    Links,
    ForceGraphIntegration,
)


class TestNodes(unittest.TestCase):
    def test_create(self):
        m_node = Nodes(id="Test_metabolite", group="metabolite")
        self.assertIsInstance(m_node, Nodes)
        self.assertEqual("Test_metabolite", m_node.id)
        self.assertEqual("metabolite", m_node.group)

        r_node = Nodes(id="Test_reaction", group="reaction")
        self.assertIsInstance(r_node, Nodes)
        self.assertEqual("Test_reaction", r_node.id)
        self.assertEqual("reaction", r_node.group)

    def test_to_json(self):
        m_node = Nodes(id="Test_metabolite", group="metabolite")
        r_node = Nodes(id="Test_reaction", group="reaction")

        self.assertEqual(
            '{"id":"Test_metabolite", "group":"metabolite"}', m_node.to_json()
        )
        self.assertEqual(
            '{"id":"Test_reaction", "group":"reaction"}', r_node.to_json()
        )


class TestLinks(unittest.TestCase):
    def test_create(self):
        link = Links(source="sourceMet", target="targetMet", value=4)
        self.assertIsInstance(link, Links)
        self.assertEqual("sourceMet", link.source)
        self.assertEqual("targetMet", link.target)
        self.assertEqual(4, link.value)

        with self.assertRaises(AttributeError):
            link.source = "anotherSource"

        with self.assertRaises(AttributeError):
            link.target = "anotherTarget"

        link.value = 7
        self.assertEqual(7, link.value)

    def test_to_json(self):
        link = Links(source="sourceMet", target="targetMet", value=4)
        self.assertEqual(
            '{"source":"sourceMet","target":"targetMet","value":4}',
            link.to_json(),
        )

        link.value = 8
        self.assertEqual(
            '{"source":"sourceMet","target":"targetMet","value":8}',
            link.to_json(),
        )

    def test_compare(self):
        link_1 = Links(source="sourceMet", target="targetMet", value=4)
        link_2 = Links(source="sourceMet", target="targetMet", value=4)
        self.assertEqual(link_1, link_2)

        link_1.value = 9
        self.assertNotEqual(link_1, link_2)

        link_2 = Links(source="sourceMet", target="targetMet", value=9)
        self.assertEqual(link_1, link_2)

        link_1 = Links(source="sourceMet", target="targetMet", value=4)
        link_2 = Links(source="sourceMet", target="targetMet2", value=4)
        self.assertNotEqual(link_1, link_2)

        link_1 = Links(source="sourceMet", target="targetMet", value=4)
        link_2 = Links(source="sourceMet2", target="targetMet", value=4)
        self.assertNotEqual(link_1, link_2)

    def test_hash(self):
        link_1 = Links(
            source="sourceMet", target="targetMet", value=4
        ).__hash__()
        link_2 = Links(
            source="sourceMet", target="targetMet", value=4
        ).__hash__()
        self.assertEqual(link_1, link_2)

        link_2 = Links(
            source="sourceMet", target="targetMet", value=8
        ).__hash__()
        self.assertEqual(link_1, link_2)

        link_2 = Links(
            source="sourceMet", target="targetMet2", value=8
        ).__hash__()
        self.assertNotEqual(link_1, link_2)

        link_2 = Links(
            source="sourceMet2", target="targetMet", value=8
        ).__hash__()
        self.assertNotEqual(link_1, link_2)

        link_2 = Links(
            source="sourceMet2", target="targetMet2", value=8
        ).__hash__()
        self.assertNotEqual(link_1, link_2)


class TestGraphData(unittest.TestCase):
    def test_create(self):
        data = GraphData()
        self.assertIsInstance(data, GraphData)
        self.assertIsInstance(data.nodes, set)
        self.assertIsInstance(data.links, set)


class TestEscher(unittest.TestCase):
    def test_create(self):
        f_graph = ForceGraphIntegration()
        self.assertIsInstance(f_graph, ForceGraphIntegration)

    def test_model(self):
        f_graph = ForceGraphIntegration()
        test_model = textbook_biocyc.copy()
        reac1 = test_model.reactions[0]
        reac2 = test_model.reactions[1]

        test_model.groups.add(
            cobra.core.Group(id="testGroup", members=[reac1, reac2])
        )

        group = test_model.groups.get_by_id("testGroup")

        f_graph.model = group
        self.assertEqual(group, f_graph.model)

    def test_solution(self):
        f_graph = ForceGraphIntegration()
        test_model = textbook_biocyc.copy()
        reac1 = test_model.reactions[0]
        reac2 = test_model.reactions[1]

        test_model.groups.add(
            cobra.core.Group(id="testGroup", members=[reac1, reac2])
        )

        solution = {"ACALD": 7}
        self.assertEqual(f_graph.solution, None)

        f_graph.solution = solution
        self.assertEqual(f_graph.solution, solution)

    def test_model_rep(self):
        f_graph = ForceGraphIntegration()
        test_model = textbook_biocyc.copy()
        reac1 = test_model.reactions[0]
        reac2 = test_model.reactions[1]

        test_model.groups.add(
            cobra.core.Group(id="testGroup", members=[reac1, reac2])
        )

        group = test_model.groups.get_by_id("testGroup")
        self.assertEqual(None, f_graph.model)

        f_graph.model = group

        expected = {
            "nodes": [
                {"id": "NAD_c", "group": "metabolite"},
                {"id": "ACETALD_e", "group": "metabolite"},
                {"id": "ACALD", "group": "reaction"},
                {"id": "NADH_c", "group": "metabolite"},
                {"id": "PROTON_c", "group": "metabolite"},
                {"id": "CO_A_c", "group": "metabolite"},
                {"id": "ACETYL_COA_c", "group": "metabolite"},
                {"id": "ACALDt", "group": "reaction"},
                {"id": "ACETALD_c", "group": "metabolite"},
            ],
            "links": [
                {"source": "ACALDt", "target": "ACETALD_c", "value": 1.0},
                {"source": "ACALD", "target": "NADH_c", "value": 1.0},
                {"source": "ACETALD_e", "target": "ACALDt", "value": 1.0},
                {"source": "ACETALD_c", "target": "ACALD", "value": 1.0},
                {"source": "CO_A_c", "target": "ACALD", "value": 1.0},
                {"source": "ACALD", "target": "ACETYL_COA_c", "value": 1.0},
                {"source": "ACALD", "target": "PROTON_c", "value": 1.0},
                {"source": "NAD_c", "target": "ACALD", "value": 1.0},
            ],
        }
        actual = json.loads(f_graph._model_rep)

        self.assertCountEqual(expected["nodes"], actual["nodes"])
        self.assertCountEqual(expected["links"], actual["links"])
