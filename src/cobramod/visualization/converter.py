#!/usr/bin/env python3
"""JSON creator

This module handles the convertion of strings into a proper object that can be
used to create a JSON string, which can be later used in Escher.

The main class of this module is
:func:`cobramod.visualization.converter.JsonDictionary`. This class is able to
parse and store data as JSON objects. To check the attributes for each
JSON object, please read the documentation of
:func:`cobramod.visualization.items`

Important methods:

- json_dump: The data can be parsed into a JSON. (str)
- add_metabolite: Add metabolite-node into the JsonDictionary.
- add_marker: Add a marker-node into the JsonDictionary.
- create_reaction: Returns a :func:`cobramod.visualization.items.Reaction`.
- add_reaction: Parses a reaction string and add the information into the
JsonDictionary.
- visualize: Saves Escher visualization as a HTML and return the Escher
Builder.
"""

from collections import UserDict
from contextlib import suppress
from itertools import cycle, chain, repeat
from json import dumps
from typing import Dict, List
from pathlib import Path
from warnings import warn, catch_warnings, simplefilter
from webbrowser import open as web_open

from escher import Builder
from IPython.core.getipython import get_ipython

from cobramod.visualization.pair import PairDictionary
from cobramod.visualization.items import Reaction, Node
from cobramod.visualization.debug import debug_log


def _in_notebook() -> bool:
    """
    Returns true if code is being executed through the IPython kernel ZMQ.
    """
    try:
        return get_ipython().__class__.__name__ == "ZMQInteractiveShell"
    except NameError:
        return False


def _convert_string(string: str) -> dict:
    """
    Converts a :func:`cobra.Reaction` reaction-string into a dictionary, and
    returns a dictionary with the corresponding participants and their
    coefficients.

    Examples:
    'C00001_c + 2 C00002_c --> C00009_c + C00080_c + G11113_c'
    'C00002_c + C00033_c <=> C00227_c + G11113_c'
    'C00002_c + C00033_c <-- C00227_c + G11113_c'
    """
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
        # FACTOR defines product or reactant
        FACTOR = 1
        if side is left:
            FACTOR = -1
        # Sometimes, there are no metabolites (sink reactions), the sigle "1"
        # must be ignored
        with suppress(ValueError):
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
        text_labels  (dict, optional): Dictionary with the custom text inside
            the canvas.
        canvas (dict, optional): x and y position, width and height of the
            white area in Escher.

    Attributes:
        CANVAS_WIDTH (float): Width for the canvas. Defaults to 1500
        CANVAS_HEIGHT (float): Height for the canvas. Defaults to 1500
        R_WIDTH (float): Width of a reaction. Defaults to 350
        R_HEIGHT (float): Height of a reaction. Defaults to 270
        reaction_data (Dict[str, float]): Dictionary with the solution to be
            visualized. Default to empty dictionary.
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
            # Check if key can be called individually
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
        self.CANVAS_WIDTH: float = self.data["canvas"]["width"]
        self.CANVAS_HEIGHT: float = self.data["canvas"]["height"]
        self.R_WIDTH: float = 350
        self.R_HEIGHT: float = 300  # 210
        # Data stored about reactions and participants.
        self._overview = dict()
        # Default solution
        self.reaction_data: Dict[str, float] = None

    def json_dump(self, indent: int = None) -> str:
        """
        Returns a string that is the JSON representation of this class.

        Args:
            indent (int): Defines the indentation for the JSON.
                Defaults to None.
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
            # Each reaction must have its Segment changed to a native dict
            for reaction in reactions:
                # if empty, must be changed to regular dictionary
                if not reactions[reaction]["segments"]:
                    reactions[reaction]["segments"] = dict()
                else:
                    # Change each Segment
                    temporal_dict = dict()
                    for key, value in reactions[reaction]["segments"].items():
                        temporal_dict[key] = dict(**value)
                    reactions[reaction]["segments"] = temporal_dict

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

    def _get_set(self, item: str) -> set:
        """
        Return a set for the keys of either the reactions, the nodes or
        segments for all reactions. Options for item: "nodes", "segments",
        "reactions"
        """
        if item == "segments":
            numbers = set()
            # Get the numbers from the segments
            with suppress(KeyError):
                reactions = self.data["reactions"]
                for key in reactions.keys():
                    segments = reactions[key]["segments"]
                    numbers.update([int(index) for index in segments])
        elif item in ("reactions", "nodes"):
            numbers = {int(key) for key in self.data[item]}
        else:
            raise ValueError(
                "Argument 'item' not correct. Refer to docstring."
            )
        return numbers

    def _get_last_number(self, item: str) -> int:
        """
        Returns the largest number of the keys from either reactions, nodes, or
        segments from each reaction. Options for item: "nodes", "segments",
        "reactions"
        """
        # Return 0 for first item, otherwise the longest number + 1
        numbers = self._get_set(item=item)
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
        Add a metabolite-type node into the JsonDictionary. The key will be
        always the last node number.

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
        number = str(self._get_last_number(item="nodes"))
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
        debug_log.info(
            f'New metabolite-node for "{bigg_id}" with id "{number}" added'
            f" to the JsonDictionary."
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
        number = str(self._get_last_number(item="nodes"))
        self.data["nodes"][number] = Node(node_type=node_type, x=x, y=y)
        debug_log.info(
            f'New {node_type}-node with id "{number}" added to '
            f"the JsonDictionary"
        )

    def __get_col_row(self) -> tuple:
        """
        Returns the number of columns and rows, in that order
        """
        columns = int(self.CANVAS_WIDTH / (self.R_WIDTH + 100))
        rows = int(self.CANVAS_HEIGHT / (self.R_HEIGHT))
        return columns, rows

    def _get_edges(self) -> tuple:
        """
        Return the top edge and the left edge for the Reaction-box. They define
        the position of the Reaction-box and takes into consideration the
        number of reactions in the JsonDictionary.
        """
        # Size of columns depends on number of reactions, or the relationship
        # CANVAS_WIDTH:R_WIDTH
        # The 50 and 100 px is a visual help (extra separation)
        # TODO: verify visual help
        columns, rows = self.__get_col_row()
        # create Generators for the place number of the reaction-Box
        sequence_x = cycle(range(0, columns))
        repeated_rows = chain.from_iterable(
            (
                repeat(row_number, columns)
                for row_number in range(0, len(self.data["reactions"]))
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
        segments: dict,
        gene_reaction_rule: str = "",
        genes: List[Dict[str, str]] = [],
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
        # TODO: Confirm behaviour of PairDictionaries
        # if not isinstance(segments, PairDictionary):
        #     raise TypeError("Argument 'segments' must be a PairDictionary")
        top_edge, left_edge = self._get_edges()
        reaction = Reaction(
            name=name,
            bigg_id=identifier,
            reversibility=reversibility,
            # label_y=top_edge + self.R_HEIGHT / 2 - 30,
            label_y=top_edge + 30,
            # label_x=left_edge + self.R_WIDTH / 2 + 20,
            label_x=(left_edge + self.R_WIDTH / 2) - len(identifier) / 2 * 17,
            gene_reaction_rule=gene_reaction_rule,
            genes=genes,
            segments=segments,
        )
        return reaction

    def _add_reaction_markers(self, identifier: str):
        """
        Add the corresponding midmarker and multimarkers into the
        JsonDictionary. Nodes will be added to the corresponding reaction data.
        Number for the nodes will not repeat themselves.
        """
        # For markers: 20 px between each one. Sequence should follow:
        # multimarker-midmarker-multimarker
        top_edge, left_edge = self._get_edges()
        for node_type, extra_x, position in (
            ("multimarker", -20, "first"),
            ("midmarker", 0, "middle"),
            ("multimarker", 20, "last"),
        ):
            last = self._get_last_number(item="nodes")
            self._overview[identifier][position] = last
            self.add_marker(
                x=left_edge + (self.R_WIDTH / 2) + extra_x + 30,
                y=top_edge + (self.R_HEIGHT / 2),
                node_type=node_type,
            )

    def _first_column(self, x_position: float) -> bool:
        """
        Returns true if given x-position is located in the first column
        """
        columns, _ = self.__get_col_row()
        x_range = self.CANVAS_WIDTH / columns
        return x_position < x_range

    def _last_column(self, x_position: float) -> bool:
        """
        Returns true if given x-position is located in the last column
        """
        columns, _ = self.__get_col_row()
        x_range = self.CANVAS_WIDTH / columns
        edge_range = self.CANVAS_WIDTH - x_range
        return x_position >= edge_range

    def _not_edges(self, x_position: float) -> bool:
        """
        Returns True if x_position is located in the last column and previous
        metabolite from the JsonDictionary is located in the first column.
        """
        # First condition/ actual reaction
        first = self._first_column(x_position=x_position)
        # Second condition/ prior reaction. Grab last metabolite, which is 4
        # identifiers behind
        last = str(self._get_last_number(item="nodes") - 4)
        second = self._last_column(x_position=self.data["nodes"][last]["x"])
        return all([first, second])

    def _add_metabolites(self, metabolite_dict: dict, reaction: Reaction):
        """
        Adds the metabolites from the dictionary into a
        :func:`cobramod.visualization.items.Reaction` and creates Nodes of the
        metabolites for the JsonDictionary class.
        """
        top_edge, left_edge = self._get_edges()
        # Minimum number of identifiers. TODO: verify behaviour with 0
        side_dict = {"left": 1, "right": 1}
        # Define previous reaction
        try:
            previous = str(len(self.data["reactions"]) - 1)
            old_name = self.data["reactions"][previous]["bigg_id"]
            old_metabolites = [
                reaction["bigg_id"]
                for reaction in self.data["reactions"][previous]["metabolites"]
                if reaction["coefficient"] > 0
            ]
            not_edges = self._not_edges(x_position=left_edge)
        except KeyError:
            # Special case for first reaction
            previous = "0"
            old_name = ""
            old_metabolites = []
            not_edges = False
        for key, value in metabolite_dict.items():
            # By default, left side
            item = "reactants"
            SIDE = 0
            counter = side_dict["left"]
            number_metabolites = len(
                [value for value in metabolite_dict.values() if value < 0]
            )
            if value > 0:
                item = "products"
                SIDE = 1
                counter = side_dict["right"]
                number_metabolites = len(
                    [value for value in metabolite_dict.values() if value > 0]
                )
            # Defining positions for metabolites based on Reaction-box
            space_y = self.R_HEIGHT / (number_metabolites + 1)
            # For Reaction-Box
            dot_x = left_edge + self.R_WIDTH * SIDE + 30
            dot_y = top_edge + counter * space_y
            # For metabolites: the labels must be 40 px higher in th y-axis
            # and the x-axis varies depending in the length (between 31-50)
            _up_or_down = cycle([-1, 1])
            if number_metabolites > 2:
                _up_or_down = cycle([-1])
            # Check previous reaction (must be in product). Side must be left
            # (reactants) AND should not be in the last column. Change last
            # variable to node of the node.
            if key in old_metabolites and SIDE == 0 and not not_edges:
                last = self._overview[old_name]["nodes"][key]
                debug_log.debug(
                    f'Metabolite "{key}" in previous reaction "{old_name}" '
                    f'located in node "{last}"'
                )
                # TODO: Move previous node
            else:
                # Add node metabolite-node to JsonDictionary. This has to be
                # omitted if previous reaction has a shared metabolite. Also,
                # change last variable to last node.
                self.add_metabolite(
                    x=dot_x,
                    y=dot_y,
                    label_x=dot_x - 30,
                    label_y=dot_y + 30 * next(_up_or_down) + 10,
                    bigg_id=key,
                    # TODO: Add proper name
                    name=key,
                    node_is_primary=False,
                )
                last = self._get_last_number(item="nodes") - 1
            # Increase number. TODO: find a better way
            counter += 1
            if value < 0:
                side_dict["left"] = counter
            else:
                side_dict["right"] = counter
            # Add coefficient to reaction object
            reaction.add_metabolite(bigg_id=key, coefficient=value)
            # Add to node dictionary. Last minus one, since the node was
            # already added.
            self._overview[reaction["bigg_id"]][item].update({key: value})
            self._overview[reaction["bigg_id"]]["nodes"].update({key: last})

    def _add_segments(self, reaction: Reaction, metabolite_dict: dict):
        # First 2 Segmenst joins the node-markers. The number of Segments is
        # equal to: 2 + number_metabolites
        identifier = reaction["bigg_id"]
        last = self._get_last_number(item="segments") - 1
        marker = {
            "first": self._overview[identifier]["first"],
            "last": self._overview[identifier]["last"],
            "middle": self._overview[identifier]["middle"],
        }
        # From markers. They will be always 2.
        for node in ("first", "last"):
            last += 1
            reaction.add_segment(
                identifier=str(last),
                from_node_id=str(marker[node]),
                to_node_id=str(marker["middle"]),
            )
        for key, value in metabolite_dict.items():
            # Two due to the first two segments, and plus one as it represent
            # the actual Segment
            number = self._overview[identifier]["nodes"][key]
            last += 1
            if value < 0:
                reaction.add_segment(
                    identifier=str(last),
                    from_node_id=str(number),
                    to_node_id=str(marker["first"]),
                )
            elif value > 0:
                reaction.add_segment(
                    identifier=str(last),
                    from_node_id=str(number),
                    to_node_id=str(marker["last"]),
                )
        # Verify the number of Segments. They cannot have the same key of other
        # Segments from other reactions.

    def __check_out_canvas(self) -> bool:
        """
        Returns true if next reaction would be within range of the canvas
        """
        rows, columns = self.__get_col_row()
        maximum = rows * columns
        actual = len(self.data["reactions"]) + 1
        return actual > maximum

    def add_reaction(self, string: str, identifier: str):
        """
        Parses and add given reaction string as a reaction for the
        JsonDictionary. It will automatically create all the necessary nodes
        and segments for the JSON.

        Args:
            string (str): Reaction string to be parsed.
            identifier (str): Identifier for the reaction

        Raises:
            UserWarning: If reaction would be located outside the canvas. It
                will not stop the method.
        """
        # Check for reaction inside canvas.
        if self.__check_out_canvas():
            msg = f'Reaction "{identifier}" will be located ouside the canvas.'
            warn(message=msg, category=UserWarning)
            debug_log.warning(msg=msg)
        # Add general data
        self._overview[identifier] = {
            "reactants": {},
            "products": {},
            "nodes": {},
        }
        # Extract information for new reaction, nr of metabolites (string
        # representation)
        metabolite_dict = _convert_string(string=string)
        # Base reaction
        reaction = self.create_reaction(
            # TODO: Change name
            name="test_reaction" + identifier,
            identifier=identifier,
            reversibility=True,
            segments=dict(),
        )
        # Add nodes (metabolites and markers)
        self._add_metabolites(
            metabolite_dict=metabolite_dict, reaction=reaction
        )
        self._add_reaction_markers(identifier=identifier)
        # Segments
        self._add_segments(reaction=reaction, metabolite_dict=metabolite_dict)
        # Add to JsonDictionary
        number = self._get_last_number(item="reactions")
        self.data["reactions"].update({str(number): reaction})
        debug_log.info(f'Reaction "{identifier}" added to the JsonDictionary.')

    def visualize(self, filepath: Path = None) -> Builder:
        """
        Saves the visualization of the JsonDictionary in given path as a HTML.
        Returns the builder for the JsonDictionary. If method is called in
        Jupyter or Qtconsole, it will show the embedded builder of the escher
        visualization. Else, it will open the default browser of the operating
        system and will load the previously saved HTML.

        Args:
            filepath (Path): Path for the HTML. Defaults to "pathway.html" in
                the current working directory
        Returns:
            Builder: Escher builder object for the visualization
        """
        if not filepath:
            filepath = Path.cwd().joinpath("pathway.html")
        # Create the builder. Text will make reactions only show the values
        builder = Builder(
            reaction_styles=["text"],
            map_name=self.data["head"]["map_name"],
            map_json=self.json_dump(),
        )
        # This statement is needed, otherwise, all reactions labels will
        # appear with "(nd)".
        if self.reaction_data:
            # TODO: only dict
            builder.reaction_data = self.reaction_data
        builder.save_html(filepath=filepath)
        # builder.reaction_styles = ["color"]
        debug_log.info(f'Visualization located in "{filepath}"')
        # If in Jupyter, launch embedded widget. Otherwise, launch web-browser
        if not _in_notebook():
            # The context manager removes the ResourceWarning
            with catch_warnings():
                simplefilter(action="ignore", category=ResourceWarning)
                web_open("file://" + str(filepath))
        return builder
