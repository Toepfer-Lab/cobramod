from collections import UserDict
from contextlib import suppress
from json import dumps
from typing import Union, Dict

from cobramod.visualization.pair import PairDictionary
from cobramod.visualization.items import Reaction, Node


class JsonDictionary(UserDict):
    def __init__(
        self,
        head: Dict[str, str] = None,
        reactions: Dict[str, Reaction] = {},
        nodes: PairDictionary = PairDictionary(),
        text_labels: Dict[str, str] = {},
        canvas: Dict[str, float] = None,
    ):
        if not head:
            head = {
                "map_name": "",
                "map_id": "",
                "map_description": "",
                "homepage": "",
                "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
            }
        # Segments from reactions and nodes must share same keys.
        if not canvas:
            canvas = {"x": 0, "y": 0, "width": 1500, "height": 1500}
        kwargs = {
            key: value
            for key, value in locals().copy().items()
            if key != "self" and key != "kwargs"
        }
        super().__init__(self, **kwargs)

    def json_dump(self, indent: Union[None, int]) -> str:
        """
        Creates a string that is the representation for JSON of the class.
        It can have an indentation.
        """
        # Transform into regular dictionaries
        nodes = {
            key: dict(**value) for key, value in self.data["nodes"].items()
        }.copy()
        # In case of reactions, the objects inside them 'Segment' must be
        # transformed as well
        try:
            reactions = {
                key: dict(**value)
                for key, value in self.data["reactions"].items()
            }.copy()
            for reaction in reactions:
                for key, value in reactions[reaction]["segments"].items():
                    reactions[reaction]["segments"][key] = dict(**value)
        except KeyError:
            reactions = {}
        return dumps(
            obj=[
                self.data["head"],
                {
                    "reactions": reactions,
                    "nodes": nodes,
                    "text_labels": self.data["text_labels"],
                    "canvas": self.data["canvas"],
                },
            ],
            indent=indent,
        )

    def _get_last_number(self) -> int:
        """
        Returns the largest number of the keys from nodes, and segments from
        each reaction.
        """
        numbers = {int(key) for key in self.data["nodes"]}
        with suppress(KeyError):
            reactions = self.data["reactions"]
            for key in reactions.keys():
                segments = reactions[key]["segments"]
                numbers.update([int(index) for index in segments])
        if not numbers:
            return 0
        return max(numbers) + 1

    def add_metabolite(
        self,
        x: float,
        y: float,
        label_x: float,
        label_y: float,
        bigg_id: str,
        name: str,
        node_is_primary: bool,
    ):
        number = str(self._get_last_number())
        self.data["nodes"][number] = Node(
            node_type="metabolite",
            x=x,
            y=y,
            bigg_id=bigg_id,
            label_x=label_x,
            label_y=label_y,
            name=name,
            node_is_primary=node_is_primary,
        )

    def add_marker(self, x: float, y: float, node_type: str = "midmarker"):
        number = str(self._get_last_number())
        self.data["nodes"][number] = Node(node_type=node_type, x=x, y=y)
