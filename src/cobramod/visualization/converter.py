from json import dumps
from typing import Union
from cobramod.visualization.pair import PairDictionary


class JsonCobramod:
    def __init__(self):
        self.head = {
            "map_name": "",
            "map_id": "",
            "map_description": "",
            "homepage": "",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
        }
        # Reactions and node must share same keys.
        self.reactions = PairDictionary()
        self.nodes = PairDictionary(pair=self.nodes)
        self.reactions.set_pair(pair=self.reactions)
        self.text_labels = {}
        self.canvas = {"x": 0, "y": 0, "width": 1500, "height": 1500}
        self.identifier = []

    def _json_dump(self, indent: Union[None, int]) -> str:
        """
        Creates a string that is the representation for JSON of the class.
        It can have an indentation.
        """
        return dumps(
            obj=[
                self.head,
                {
                    "reactions": self.reactions,
                    "nodes": self.nodes,
                    "text_labels": self.text_labels,
                    "canvas": self.canvas,
                },
            ],
            indent=indent,
        )

    def _check_identifier(self):
        pass
