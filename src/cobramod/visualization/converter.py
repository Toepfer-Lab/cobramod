from collections import UserDict
from contextlib import suppress
from itertools import cycle, chain, repeat
from json import dumps
from typing import Dict, List

from cobramod.visualization.pair import PairDictionary
from cobramod.visualization.items import Reaction, Node


def _convert_string(string: str) -> dict:
    # C00001_c + 2 C00002_c --> C00009_c + C00080_c + G11113_c
    # 'C00002_c + C00033_c <=> C00227_c + G11113_c'
    # 'C00002_c + C00033_c <-- C00227_c + G11113_c'
    middle = max(string.find(">"), string.find("<"))
    # find exact middle
    if " " == string[middle - 1]:
        middle = middle + 1
    else:
        middle = middle - 1
    left, right = string[: middle - 1], string[middle + 2 :]
    metabolites = dict()
    for side in (left, right):
        # If split do not work, then add a "1"
        sequence = [
            item.split() if len(item.split()) == 2 else ["1", *item.split()]
            for item in side.split("+")
        ]
        FACTOR = 1
        if side is left:
            FACTOR = -1
        metabolites.update(
            {key: float(value) * FACTOR for value, key in sequence}
        )
    return metabolites


class JsonDictionary(UserDict):
    """
    Create a JsonDictionary object which can be used to parse information into
    JSON files to be later read by Escher. When creating the object, some
    keyword arguments may be passed.

    Keyword Arguments:
        head (dict, optional): General information of the JSOB. Keys
            included: map_name, map_id, map_description, homepage, schema
        reactions (dict, optional): Dictionary with multiple
            :func:`cobramod.visualization.items.Reaction` where the key is the
            number of the reaction and the value the object. Defaults to empty
            dictionary
        nodes (PairDictionary): Dictionary with multiple
            :func:`cobramod.visualization.items.Node` where the key is the
            number of the Node and the value the corresponding object. Defaults
            to empty dictionary.
        text_labels  (dict, optional): Dictionary with the custom text in the
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
        # Defining variables for sizes
        self.CANVAS_WIDTH = self.data["canvas"]["width"]
        self.CANVAS_HEIGHT = self.data["canvas"]["height"]
        self.R_WIDTH = 350
        self.R_HEIGHT = 210

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
                # if empty, must be changed to regular dictionary
                if not reactions[reaction]["segments"]:
                    reactions[reaction]["segments"] = dict()
                else:
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

    def create_marker(self, x: float, y: float, node_type: str = "midmarker"):
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

    def _get_edges(self) -> tuple:
        """
        Return the top edge and the left edge for the Reaction-box. They define
        the position of the Reaction-box and takes into consideration the
        number of reactions in the JsonDictionary.
        """
        # Size of columns depends on number of reactions, or the relationship
        # CANVAS_WIDTH:R_WIDTH
        # The 50 and 80 px is a visual help (extra separation)
        columns = min(
            [
                # len(self.data["reactions"]) + 2,
                int(self.CANVAS_WIDTH / (self.R_WIDTH + 50))
            ]
        )
        rows = min(
            [
                # ADD: minium value from columns,
                int(self.CANVAS_HEIGHT / (self.R_HEIGHT + 80))
            ]
        )
        # create Generators for the place number of the reaction-Box
        sequence_x = cycle(range(0, columns))
        repeated_rows = list(
            chain.from_iterable(
                (
                    repeat(row_number, columns)
                    for row_number in range(0, len(self.data["reactions"]))
                )
            )
            # A 0 just in case of first item
        ) or [0]
        sequence_y = cycle(repeated_rows)
        length = list(range(0, len(self.data["reactions"])))
        # TODO: find more elegant way
        if not length:
            place_x = 0
            place_y = 0
        else:
            # Skip once 0 for both cases because of if-statement
            next(sequence_x)
            next(sequence_y)
            for repetition in length:
                place_x = next(sequence_x)
                place_y = next(sequence_y)
        # External box
        top_edge = (self.CANVAS_HEIGHT / rows) * place_y
        left_edge = (self.CANVAS_WIDTH / columns) * place_x
        return top_edge, left_edge

    def create_reaction(
        self,
        name: str,
        identifier: str,
        reversibility: bool,
        gene_reaction_rule: str = "",
        genes: List[Dict[str, str]] = [],
        segments: PairDictionary = PairDictionary(),
    ) -> Reaction:
        """
        Returns a :func:`cobramod.visualization.items.Reaction`. It will take
        into consideration the actual number of reactions in the JsonDictionary
        and add proper x and y position for the labels.

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
        top_edge, left_edge = self._get_edges()
        reaction = Reaction(
            name=name,
            bigg_id=identifier,
            reversibility=reversibility,
            label_y=top_edge + self.R_HEIGHT / 2 - 30,
            label_x=left_edge + self.R_WIDTH / 2,
            gene_reaction_rule=gene_reaction_rule,
            genes=genes,
            segments=segments,
        )
        return reaction

    def _add_metabolites(self, metabolite_dict: dict, reaction: Reaction):
        """
        Adds the metabolites from the dictionary into a
        :func:`cobramod.visualization.items.Reaction` and creates Nodes of the
        metabolites for the JsonDictionary class.
        """
        top_edge, left_edge = self._get_edges()
        # Minimum number of identifiers. TODO: verify behaviour with 0
        side_dict = {"left": 1, "right": 1}
        # Internal box (Reaction-box)
        for key, value in metabolite_dict.items():
            # By default, right side
            SIDE = 0
            counter = side_dict["right"]
            number_metabolites = len(
                [value for value in metabolite_dict.values() if value > 0]
            )
            if value < 0:
                SIDE = 1
                counter = side_dict["left"]
                number_metabolites = len(
                    [value for value in metabolite_dict.values() if value < 0]
                )
            # Defining positions for metabolites based on Reaction-box
            space_y = self.R_HEIGHT / (number_metabolites + 1)
            # For Reaction-Box
            dot_x = left_edge + self.R_WIDTH * SIDE + 30
            dot_y = top_edge + counter * space_y
            # For metabolites: the labels must be 40 px higher in th y-axis
            # and the x-axis varies depending in the length (between 31-50)
            self.add_metabolite(
                x=dot_x,
                y=dot_y,
                label_x=dot_x - 30,
                label_y=dot_y + 40,
                bigg_id=key,
                name="test_" + key,
                node_is_primary=False,
            )
            # Increase number. TODO: find a better way
            counter += 1
            if value < 0:
                side_dict["left"] = counter
            else:
                side_dict["right"] = counter
            reaction.add_metabolite(bigg_id="test_" + key, coefficient=value)

    def _add_markers(self):
        # For markers: 20 px between each one. Sequence should follow:
        # multimarker-midmarker-multimarker
        top_edge, left_edge = self._get_edges()
        for node_type, extra_x in (
            ("multimarker", -20),
            ("midmarker", 0),
            ("multimarker", 20),
        ):
            self.create_marker(
                x=left_edge + (self.R_WIDTH / 2) + extra_x + 30,
                y=top_edge + (self.R_HEIGHT / 2),
                node_type=node_type,
            )

    def shift(self):
        row_size = int(self.CANVAS_WIDTH / (self.R_WIDTH + 50))
        # if row_size > actual_length, shift only reactions that surpassed it.
        if row_size > len(self.data["reactions"]):
            for reaction in self.data["reactions"].values():
                reaction["label_x"] += -400
                print(".")

    def add_reaction(self, string: str, identifier: str):
        # Extract information for new reaction, nr of metabolites (string
        # representation)
        metabolite_dict = _convert_string(string=string)
        # if max_size is not surpassed, nothing changes
        # shift if needed
        if len(self.data["reactions"]) >= 1:
            # self.shift()
            # Modify prior reactions (position)
            pass
        # Base reaction
        reaction = self.create_reaction(
            name="test_reaction" + identifier,
            identifier="test_R" + identifier,
            reversibility=True,
        )
        # Add nodes
        self._add_metabolites(
            metabolite_dict=metabolite_dict, reaction=reaction
        )
        self._add_markers()
        # Segments
        # Verify the number of Segments. They cannot have the same key of other
        # Segments from other reactions.
        # Add to JsonDictionary
        number = self._get_last_number(reaction=True)
        self.data["reactions"].update({str(number): reaction})
