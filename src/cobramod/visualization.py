from contextlib import suppress
from json import dumps
from typing import Union, Any

from cobramod.error import NodeAttributeError


class Node:
    def __new__(
        self,
        node_type: str,
        x: float,
        y: float,
        label_x: float = None,
        label_y: float = None,
        bigg_id: str = "",
        name: str = "",
        node_is_primary: bool = False,
    ) -> Any:
        """
        Build a node for JsonCobramod, which can be parsed later into a proper
        JSON for Escher. A node can be a metabolite, midmarker or multimarker.
        Markers are the visible dots that are located in the middle of the
        reaction that defines the orientation (vertical, horizonzal).
        Metabolites represent the orange dots.


        Args:
            node_type: Can either a 'metabolite', 'midmarker' or 'multimarker'.
                Metabolites must have all arguments. Markers only need y and x.
            x (float): Location in x-axis for the node.
            y (float): Location in y-axis for the node
            label_x (float): Location in x-axis for the label of the node.
            label_y (float): Location in x-axis for the label of the node.
            bigg_id (str, optional): Identifier for the metabolite. It does not
                have to be a bigg identifier.
            name (str, optional): The name for the metabolite.
            node_is_primary (bool, optional): If the metabolite is a primary
                compound.

        Returns:
            dict: Dictionary with the information of the type of node.
        """
        # Using a copy of locals, without self and the kwargs
        kwargs = {
            key: value
            for key, value in locals().copy().items()
            if key != "self" and key != "kwargs"
        }
        for method in (
            self.metabolite(**kwargs),
            self.marker(node_type=node_type, x=x, y=y),
        ):
            with suppress(NodeAttributeError, TypeError):
                return next(method)
        raise NodeAttributeError(
            "Wrong attributes. Check 'node_types' and fix the arguments."
        )

    @staticmethod
    def marker(node_type: str, x: float, y: float):
        """
        Returns a dictionary with the arguments if node_type correspond to
        midmarker or multimarker.
        """
        # Must be yield, otherwise, the for loop wont work
        if node_type in ("midmarker", "multimarker"):
            yield locals().copy()
        raise NodeAttributeError

    @staticmethod
    def metabolite(
        node_type: str,
        x: float,
        y: float,
        label_x: float,
        label_y: float,
        bigg_id: str = "",
        name: str = "",
        node_is_primary: bool = False,
    ):
        """
        Returns a dictionary with the arguments if node_type correspond to
        metabolite.
        """
        # Must be yield, otherwise, the for loop wont work
        if node_type == "metabolite":
            yield locals().copy()
        raise NodeAttributeError


class JsonCobramod:
    def __init__(self):
        self.head = {
            "map_name": None,
            "map_id": None,
            "map_description": None,
            "homepage": None,
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
        }
        self.reactions = {}
        self.nodes = {}
        self.text_labels = {}
        self.canvas = {"x": 0, "y": 0, "width": 1500, "height": 1500}

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


json_test = JsonCobramod()
A = Node(node_type="metabolite", x=150, y=150, label_x=1, label_y=2)
B = Node(node_type="midmarker", x=150, y=150)
json_test.nodes["0"] = A
json_test.nodes["1"] = B
print(json_test._json_dump(indent=4))
