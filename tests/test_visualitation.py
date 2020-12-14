from pathlib import Path
import unittest
from logging import DEBUG

from cobramod.debug import debug_log
from cobramod.error import FoundInPairError
from cobramod.visualization.pair import PairDictionary
from cobramod.visualization.items import Node, Segment, Reaction
from cobramod.visualization.converter import JsonDictionary

debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")


class TestVisualization(unittest.TestCase):
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

    def test_JsonDictionary(self):
        # CASE 0: Checking behaviour with two instances
        test_dict = JsonDictionary()
        test_dict_2 = JsonDictionary()
        test_dict["reactions"]["0"] = "test_string"
        self.assertIsNot(expr1=test_dict, expr2=test_dict_2)
        self.assertRaises(KeyError, lambda: test_dict_2["reactions"]["0"])
        # CASE 1a: creation of dictionary without extra arguments.
        test_dict = JsonDictionary()
        self.assertEqual(first={}, second=test_dict["reactions"])
        self.assertEqual(first="", second=test_dict["head"]["map_name"])
        self.assertEqual(first=1500, second=test_dict["canvas"]["width"])
        # CASE 2: get last number in JsonDictionary. Reactions are not included
        # but their segments
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
            first=2, second=test_dict._get_last_number(reaction=False)
        )
        # CASE 3a: method create_reaction. Regular.
        test_dict = JsonDictionary()
        test_dict.create_reaction(
            name="test_reaction", bigg_id="test_id", reversibility=True
        )
        self.assertEqual(
            first="test_id", second=test_dict["reactions"]["0"]["bigg_id"]
        )
        # CASE 3b: checking index in dictionaries
        test_dict.create_reaction(
            name="test_reaction_2", bigg_id="test_id_2", reversibility=True
        )
        self.assertEqual(
            first="test_id_2", second=test_dict["reactions"]["1"]["bigg_id"]
        )

    def test_automate(self):
        test_dict = JsonDictionary()
        test_dict.add_reaction()
        with open(file="test.json", mode="w+") as f:
            f.write(test_dict.json_dump(indent=4))
        print(test_dict.json_dump(indent=4))
        pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
