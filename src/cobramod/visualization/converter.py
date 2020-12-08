from collections import UserDict
from contextlib import suppress
from json import dumps
from typing import Dict, List

from cobramod.visualization.pair import PairDictionary
from cobramod.visualization.items import Reaction, Node


class JsonDictionary(UserDict):
    """
    Create a JsonDictionary object which can be used to parse information into
    JSON files to be later read by Escher. When creating the object, some
    keyword arguments may be passed.

    Keyword Arguments:
        head: dict[str, str] = none,
        reactions: dict[str, reaction] = {},
        nodes: pairdictionary = pairdictionary(),
        text_labels: dict[str, str] = {},
        canvas: dict[str, float] = none,
    """

    def __init__(self, *args, **kwargs):
        """
        Initiates the creation of the JsonDictionary. It uses the same
        arguments and keyword arguments as a regular dictionary. However, some
        important keys can be imported, shown in the docstring for the class.
        """
        # Init dictionary
        super().__init__(self, *args, **kwargs)
        # Check if no kwargs were specified
        for key in ("head", "reactions", "text_labels", "nodes", "canvas"):
            # Check if key can be called individualley
            try:
                self.data[key]
            except KeyError:
                # If arguments are expected in the method, initialization of
                # objects gets bugged. e.g. modifying attribute in class A,
                # would change the same attribute in Class B
                if key == "head":
                    self.data[key] = {
                        "map_name": "",
                        "map_id": "",
                        "map_description": "",
                        "homepage": "",
                        "schema": "https://escher.github.io/escher/"
                        + "jsonschema/1-0-0#",
                    }
                elif key == "reactions":
                    self.data[key] = {}
                elif key == "nodes":
                    self.data[key] = PairDictionary()
                elif key == "text_labels":
                    self.data[key] = {}
                elif key == "canvas":
                    self.data[key] = {
                        "x": 0,
                        "y": 0,
                        "width": 1500,
                        "height": 1500,
                    }

    def json_dump(self, indent: int = None) -> str:
        """
        Creates a string that is the representation for JSON of the class.
        It can have indentation.
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

    def _get_set(self, reaction: bool):
        """
        Returns set for either the keys of the reactions or, the keys from
        nodes and the corresponding segements from the reactions. Argument
        'reaction' will return the index of reactions if true.
        """
        if not reaction:
            numbers = {int(key) for key in self.data["nodes"]}
            # Get the numbers from the segments
            with suppress(KeyError):
                reactions = self.data["reactions"]
                for key in reactions.keys():
                    segments = reactions[key]["segments"]
                    numbers.update([int(index) for index in segments])
        else:
            numbers = {int(key) for key in self.data["reactions"]}
        return numbers

    def _get_last_number(self, reaction: bool) -> int:
        """
        Returns the largest number of the keys from nodes, and segments from
        each reaction.
        """
        # Return 0 for first item, otherwise the longest number + 1
        numbers = self._get_set(reaction=reaction)
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
        """
        Add a metabolite-type node into the JsonDictionary.

        Args:
            x (float): Position in x-axis for the node.
            y (float): Position in y-axis fot the node.
            label_x (float): Position of the label in the x-axis.
            label_y (float): Position of the label in the y-axis.
            bigg_id (str): Identifier of the metabolite. It does not have to be
                from BIGG.
            name (str): Name of the label for the metabolite.
            node_is_primary (bool): True if node should represent a primary
                metabolite, i.e. Node is visually larger.

        """
        number = str(self._get_last_number(reaction=False))
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
        """
        Add a marker-type node into the JsonDictionary. Node can be a midmarker
        or a multimarker. These markes are located in the middle of the
        reaction. A midmarker is located between two multimarkers.

        Args:
            x (float): Position in x-axis for the node.
            y (float): Position in y-axis for the node.
            node_type (str): Type of marker. Options: 'midmarker' or
                'multimarker'
        """
        number = str(self._get_last_number(reaction=False))
        self.data["nodes"][number] = Node(node_type=node_type, x=x, y=y)

    def add_reaction(
        self,
        name: str,
        bigg_id: str,
        reversibility: bool,
        label_y: float = 0,
        label_x: float = 0,
        gene_reaction_rule: str = "",
        genes: List[Dict[str, str]] = [],
        segments: PairDictionary = PairDictionary(),
    ):
        """
        Adds a :func:`cobramod.visualization.items.Reaction` into the
        JsonDictionary. The segments of the reaction will be paired to the
        nodes of the JsonDictionary.

        Args:
            name (str): The name for the reaction.
            bigg_id (str): Identifier for the reaction. It does not have to be
                from BIGG.
            reversibility (bool): True if the reaction should be represented as
                reversible.
            label_y (float, optional): Location in x-axis for the label of the
                reaction. Defaults to 0.
            label_x (float, optional): Location in y-axis for the label of the
                reaction. Defaults to 0.
            gene_reaction_rule (str, optional): Gene rules that specify
                involved genes. Defaults to an empty string.
            genes (list, optional): A list with the genes involved for that
                identifier. Each item has to be a dictionary with the keys
                'bigg_id' and 'name'
            segments (PairDictionary, optional): Dictionary with segments,
                which represent the conections between nodes.
        """
        if not isinstance(segments, PairDictionary):
            raise TypeError("Argument 'segments' must be a PairDictionary")
        reaction = Reaction(
            name=name,
            bigg_id=bigg_id,
            reversibility=reversibility,
            label_y=label_y,
            label_x=label_x,
            gene_reaction_rule=gene_reaction_rule,
            genes=genes,
            segments=segments,
        )
        # Pairing and adding
        # In case of empty nodes
        reaction["segments"].set_pair(pair=self.data["nodes"])
        number = self._get_last_number(reaction=True)
        self.data["reactions"][str(number)] = reaction
