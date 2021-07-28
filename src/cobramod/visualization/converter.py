#!/usr/bin/env python3
"""JSON creator

This module handles the convertion of strings into a proper object that can be
used to create a JSON string, which can be later used in Escher.

The main class of this module is
:class:`cobramod.visualization.converter.JsonDictionary`. This class is able to
parse and store data as JSON objects. To check the attributes for each
JSON object, please read the documentation of
:mod:`cobramod.visualization.items`.

Important methods:


- json_dump: The data can be parsed into a JSON. (str)
- add_metabolite: Add metabolite-node into the JsonDictionary.
- add_marker: Add a marker-node into the JsonDictionary.
- create_reaction: Returns a :class:`cobramod.visualization.items.Reaction`.
- add_reaction: Parses a reaction string and add the information into the
JsonDictionary.
- add_blank: Adds a empty reaction. This is useful for the extra space in the
visualizations.
- visualize: Saves Escher visualization as a HTML and return the Escher
Builder.
"""
import math

import webcolors
from collections import UserDict, namedtuple
from contextlib import suppress
from itertools import cycle
from json import dumps
from typing import Dict, Optional, List, Union, Tuple
from pathlib import Path
from warnings import catch_warnings, simplefilter
from webbrowser import open as web_open

import numpy as np
from escher import Builder
from IPython.core.getipython import get_ipython

# from cobramod.visualization.pair import PairDictionary
from numpy import ndarray

from cobramod.visualization.items import Reaction, Node
from cobramod.visualization.debug import debug_log
from cobramod.visualization.mapping import get_mapping, transpose

Position = namedtuple("Position", ["row", "column"])


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
    Converts a :class:`cobra.Reaction` reaction-string into a dictionary, and
    returns a dictionary with the corresponding participants and their
    coefficients.

    Examples:
    'C00001_c + 2 C00002_c --> C00009_c + C00080_c + G11113_c'
    'C00002_c + C00033_c <=> C00227_c + G11113_c'
    'C00002_c + C00033_c <-- C00227_c + G11113_c'
    """
    # TODO: return reversibility
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
        # Sometimes, there are no metabolites (sink reactions), the single "1"
        # must be ignored
        with suppress(ValueError):
            metabolites.update(
                {key: float(value) * FACTOR for value, key in sequence}
            )
    return metabolites


def _color2np_rgb(color: Union[str, List[int], None]) -> ndarray:
    """
    This function translates the rgb or string representation into a numpy
    array. The string representation must follow the css standard.

    Args:
        color: The color to be translated.

    Returns:
        Numpy array with the rgb representation of the color or with the
        values [220,220,220] as standard return for incorrect inputs.
    """
    try:
        assert color is not None
        if not all(isinstance(elements, int) for elements in color):
            color = webcolors.name_to_rgb(color)
    except (KeyError, TypeError, ValueError, AssertionError):
        debug_log.warning(
            f'Unknown color or wrong format: "{color}". Using default color.'
        )
        color = [220, 220, 220]

    # turn int arrays into numpy arrays
    return np.array(color, dtype=np.float32)


def _divide_values(
    flux: List[float], min_max: Optional[List[float]] = None
) -> Tuple[List[float], List[float], bool, bool]:
    """
    This function divides a list into two, one consisting of the positive
    values and one of the negative values.

    Args:
        flux (List[int]): Array that is to be divided into positive and
            negative values.
        min_max ([int,int]): List consisting of two values. These values
            determine the maximum value and minimum value that are taken into
            account in the distribution. All values outside this interval
            are ignored.
    Returns:
        This function returns two lists and two bools. The first return is the
            list consisting of positive values. The second return is the list
            with negative values. The third return is a bool that indicates
            whether min_max consists of two positive values. The fourth
            return describes whether min_max consists of two negative values.
    """

    # if required, take min and max settings into account and add them if
    # not present

    both_positive = False
    both_negative = False

    for _ in [0]:
        if min_max is None:
            continue

        if min_max[0] > min_max[1]:
            debug_log.warn(
                "Set minimum is greater than maximum. Ignoring min_max"
            )
            continue

        both_positive = min_max[0] > 0 and min_max[1] > 0
        both_negative = min_max[0] < 0 and min_max[1] < 0

        size = len(flux)
        msg = (
            f"Due to set min_max values were ignored. Original range "
            f"was [{min(flux)},{max(flux)}] "
            f"set is [{min_max[0]},{min_max[1]}]."
        )

        flux = [value for value in flux if min_max[0] < value < min_max[1]]
        if size > len(flux):
            debug_log.info(msg)

        if min_max[0] != 0 and min_max[0] not in flux:
            flux.append(min_max[0])

        if min_max[1] != 0 and min_max[1] not in flux:
            flux.append(min_max[1])

    # divide positive and negative values
    positive = list()
    negative = list()

    for value in flux:
        if value > 0:
            positive.append(value)
        elif value < 0:
            negative.append(value)

    return positive, negative, both_positive, both_negative


class JsonDictionary(UserDict):
    """
    Create a JsonDictionary object which can be used to parse information into
    JSON files to be later read by Escher. When creating the object, some
    keyword arguments may be passed.

    Keyword Arguments:
        head (dict, optional): General information of the JSON. Keys
            included: map_name, map_id, map_description, homepage, schema.
        reactions (dict, optional): Dictionary with multiple
            :class:`cobramod.visualization.items.Reaction` where the key is the
            number of the reaction and the value the object. Defaults to empty
            dictionary.
        nodes (PairDictionary): Dictionary with multiple
            :class:`cobramod.visualization.items.Node` where the key is the
            number of the Node and the value the corresponding object. Defaults
            to empty dictionary.
        text_labels  (dict, optional): Dictionary with the custom text inside
            the canvas.
        canvas (dict, optional): x and y position, width and height of the
            white area in Escher.

    Attributes:
        CANVAS_WIDTH (float): Width for the canvas. Defaults to 1500.
        CANVAS_HEIGHT (float): Height for the canvas. Defaults to 1500.
        R_WIDTH (float): Width of a reaction. Defaults to 350.
        R_HEIGHT (float): Height of a reaction. Defaults to 270.
        reaction_data (Dict[str, float]): Dictionary with the solution to be
            visualized. Default to empty dictionary.
    """

    def __init__(self):
        """
        Initiates the creation of the JsonDictionary. It uses the same
        arguments and keyword arguments as a regular dictionary. However, some
        important keys can be imported, shown in the docstring for the class.
        """
        # Init dictionary
        super().__init__()
        # Initiate basic information
        self.data["head"] = {
            "map_name": "",
            "map_id": "",
            "map_description": "",
            "homepage": "",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
        }
        self.data["reactions"] = {}
        # TODO: add PairDictionary. Check if even needed
        self.data["nodes"] = {}
        self.data["text_labels"] = {}
        # Canvas variables
        self.CANVAS_WIDTH: float = 0
        self.CANVAS_HEIGHT: float = 0
        # TODO: fix variables
        self.X = 0
        self.Y = 0
        # Reaction size
        self.R_WIDTH: float = 650
        self.R_HEIGHT: float = 450  # 210
        # Data stored about reactions and participants.
        self._overview = dict()
        # Default solution
        self.flux_solution: Dict[str, float] = None
        # Dictionary with relationship of reactions
        self.graph: dict = dict()
        self.reaction_strings = dict()
        self.reaction_scale = dict()

    def get_canvas(self) -> dict:
        return {
            "x": self.X,
            "y": self.Y,
            "width": self.CANVAS_WIDTH,
            "height": self.CANVAS_HEIGHT,
        }

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
                    "canvas": self.get_canvas(),
                },
            ],
            indent=indent,
        )

    def _get_set(self, item: str) -> set:
        """
        Return a set for the keys of either the reactions, the nodes or
        segments for all reactions. Options for item: "nodes", "segments",
        "reactions".
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
        "reactions".
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

    def add_marker(self, x: float, y: float, node_type: str):
        """
        Add a marker-type node into the JsonDictionary. Node can be a midmarker
        or a multimarker. These markers are located in the middle of the
        reaction. A midmarker is located between two multimarkers.

        Args:
            x (float): Position in x-axis for the node.
            y (float): Position in y-axis for the node.
            node_type (str): Type of marker. Options: 'midmarker' or
                'multimarker'.
        """
        number = str(self._get_last_number(item="nodes"))
        self.data["nodes"][number] = Node(node_type=node_type, x=x, y=y)
        debug_log.info(
            f'New {node_type}-node with id "{number}" added to '
            f"the JsonDictionary"
        )

    def _add_reaction_markers(
        self,
        identifier: str,
        top_edge: float,
        left_edge: float,
        vertical: bool,
    ):
        """
        Add the corresponding midmarker and multimarkers into the
        JsonDictionary. Nodes will be added to the corresponding reaction data.
        Number for the nodes will not repeat themselves. If vertical is passed,
        then the reaction markers will be created in order to represent a
        vertical pathway.
        """
        # For markers: 20 px between each one. Sequence should follow:
        # multimarker-midmarker-multimarker
        # if top_edge is None and left_edge is None:
        #     top_edge, left_edge = self._get_edges()
        for node_type, extra_x, node_position in (
            ("multimarker", -20, "_first"),
            ("midmarker", 0, "_middle"),
            ("multimarker", 20, "_last"),
        ):
            last = self._get_last_number(item="nodes")
            self._overview[identifier]["nodes"][node_position] = last
            # Defining position of labels depending of position axis.
            if vertical:
                self.add_marker(
                    x=left_edge + (self.R_WIDTH / 2),
                    y=top_edge + (self.R_HEIGHT / 2) + extra_x + 0,
                    node_type=node_type,
                )
            else:
                self.add_marker(
                    x=left_edge + (self.R_WIDTH / 2) + extra_x + 30,
                    y=top_edge + (self.R_HEIGHT / 2),
                    node_type=node_type,
                )

    def _get_matrix_reactions(self, vertical: bool, position: int) -> list:
        """
        Return a List with the names of the reactions, which are located in
        given position. The position can be from a column or a row. Vertical
        defines the orientation.
        """
        if not vertical:
            data = [
                reaction
                for reaction in self._overview.keys()
                if self._overview[reaction]["position"].column == position
            ]
        else:
            data = [
                reaction
                for reaction in self._overview.keys()
                if self._overview[reaction]["position"].row == position
            ]
        return data

    def _get_products(self, reactions: list) -> Dict[str, list]:
        """
        Returns a dictionary where keys are the reactions of given list and the
        values the products of the reactions.

        Args:
            reactions (list): List with reaction identifiers.

        Returns:
            dict: Dictionary with products for given reactions.
        """
        previous = dict()
        # Find in class for the Reaction object to find their metabolites
        for reaction in self._overview.keys():
            if reaction in reactions:
                index = self._overview[reaction]["index"]
                try:
                    reaction_obj: Reaction = self.data["reactions"][index]
                    # Obtaining product side of the reactions
                    # key metabolites returns a list
                    previous[reaction] = [
                        metabolite["bigg_id"]
                        for metabolite in reaction_obj["metabolites"]
                        if metabolite["coefficient"] > 0
                    ]
                except Exception:
                    pass
        return previous

    def _find_shared(
        self, metabolite: str, products: Dict[str, list]
    ) -> tuple:
        """
        Returns the node number of the metabolite and the reaction involved if
        found in given dictionary with products.

        Args:
            metabolite (str): Identifier of the metabolite.
            products (dict): Dictionary with reactions and their corresponding
                products as values.

        Returns:
            tuple: the node_number for the metabolite and the name of the
                involved reaction.
        """
        # If nothing is found, empty strings
        node_number = str()
        old_reaction = str()
        for reaction, metabolites_list in products.items():
            if metabolite in metabolites_list:
                node_number = self._overview[reaction]["nodes"][metabolite]
                old_reaction = reaction
                # Only one is necessary
                break
        return node_number, old_reaction

    def map_metabolites(
        self,
        metabolite_dict: dict,
        reaction: Reaction,
        top_edge: float,
        left_edge: float,
        vertical: bool,
    ):
        """
        Creates the metabolites from given dictionary and complements the
        :class:`cobramod.visualization.items.Reaction`. Moreover, it creates
        the corresponding metabolites-nodes for the JsonDictionary class.

        Args:
            metabolite_dict (dict): Dictionary with metabolites and their
                coefficients.
            reaction (Reaction): Reaction that will include the metabolite.
            top_edge (float): Position for the top edge of the reaction-box.
            left_edge (float): Position for the left edge of the reaction-box.
            vertical (bool): Defines if metabolite should be map vertically.
        """
        # Obtaining position of reaction
        identifier = reaction["bigg_id"]
        position: Position = self._overview[identifier]["position"]
        # Minimum number of identifiers. TODO: verify behavior with 0
        side_dict = {"left": 1, "right": 1}
        # Obtain products from reactions in the prior column.
        # The column must be actual - 1, or 0. This is to check for shared
        # metabolites.
        position_value = position.column
        if vertical:
            position_value = position.row
        # Obtains reactions that shared either the previous column or row
        shared_reactions = self._get_matrix_reactions(
            vertical=vertical, position=max(position_value - 1, 0)
        )
        # Remove actual identifier
        filtered = [
            reaction for reaction in shared_reactions if reaction != identifier
        ]
        shared_products = self._get_products(reactions=filtered)
        # Add metabolites from dictionary.
        # TODO: get rid off if-statements
        for metabolite, coefficient in metabolite_dict.items():
            # By default, left side
            item = "reactants"
            SIDE = 0
            counter = side_dict["left"]
            number_metabolites = len(
                [value for value in metabolite_dict.values() if value < 0]
            )
            if coefficient > 0:
                item = "products"
                SIDE = 1
                counter = side_dict["right"]
                number_metabolites = len(
                    [value for value in metabolite_dict.values() if value > 0]
                )
            # For metabolites: the labels must be 40 px higher in th y-axis
            # and the x-axis varies depending in the length (between 31-50)
            _up_or_down = cycle([-1, 1])
            if number_metabolites > 2:
                _up_or_down = cycle([-1])
            # Defining positions and labels depending on the Reaction-box and
            # the axis given
            if vertical:
                space_y = self.R_WIDTH / (number_metabolites + 1)
                dot_x = left_edge + counter * space_y
                dot_y = top_edge + (self.R_HEIGHT * 0.8 * SIDE) + 30
                label_x = dot_x + 10
                label_y = dot_y - 10 + SIDE * 40
            else:
                space_y = self.R_HEIGHT / (number_metabolites + 1)
                dot_x = left_edge + (self.R_WIDTH * 0.8 * SIDE) + 30
                dot_y = top_edge + counter * space_y
                label_x = dot_x - 30
                label_y = dot_y + 30 * next(_up_or_down) + 10
            # Check for shared metabolites only if metabolite is located in
            # reactant side
            if SIDE == 0:
                shared_node, old_reaction = self._find_shared(
                    metabolite=metabolite, products=shared_products
                )
            else:
                shared_node, old_reaction = str(), str()
            # If not empty, then defined new node
            if shared_node:
                node_number = shared_node
                debug_log.debug(
                    f'Metabolite "{metabolite}" in previous reaction '
                    f'"{old_reaction}" located in node "{shared_node}"'
                )
                # TODO: Move previous node
            else:
                # Add node metabolite-node to JsonDictionary. This has to be
                # omitted if previous reaction has a shared metabolite. Also,
                # change last variable to last node.
                self.add_metabolite(
                    x=dot_x,
                    y=dot_y,
                    label_x=label_x,
                    label_y=label_y,
                    bigg_id=metabolite,
                    # TODO: Add proper name
                    name=metabolite,
                    node_is_primary=False,
                )
                node_number = self._get_last_number(item="nodes") - 1
            # Increase number. TODO: find a better way
            counter += 1
            if coefficient < 0:
                side_dict["left"] = counter
            else:
                side_dict["right"] = counter
            # Add coefficient to reaction object
            reaction.add_metabolite(
                bigg_id=metabolite, coefficient=coefficient
            )
            # Add to node dictionary. Last minus one, since the node was
            # already added.
            # Either product or reactant
            self._overview[identifier][item].update({metabolite: coefficient})
            self._overview[identifier]["nodes"].update(
                {metabolite: node_number}
            )

    def add_segments(self, reaction: Reaction, metabolite_dict: dict):
        """
        Add the segments to given Reaction. This will make the visuals in
        Escher. The information about the nodes of metabolites in located in
        the JsonDictionary. A reaction will have 2 + number of metabolite as
        its number of segments.

        Args:
            metabolite_dict (dict): Dictionary with metabolites and their
                coefficients.
            reaction (Reaction): Reaction to extend.
        """
        # Defining identifier, the last number of segments and the markers
        identifier = reaction["bigg_id"]
        last = self._get_last_number(item="segments") - 1
        marker = {
            "first": self._overview[identifier]["nodes"]["_first"],
            "last": self._overview[identifier]["nodes"]["_last"],
            "middle": self._overview[identifier]["nodes"]["_middle"],
        }
        # From markers. They will be always 2.
        for node in ("first", "last"):
            last += 1
            reaction.add_segment(
                identifier=str(last),
                from_node_id=str(marker[node]),
                to_node_id=str(marker["middle"]),
            )
        for metabolite, coefficient in metabolite_dict.items():
            # Two due to the first two segments, and plus one as it represent
            # the actual Segment
            number = self._overview[identifier]["nodes"][metabolite]
            last += 1
            # Check whether reactant or product
            if coefficient < 0:
                reaction.add_segment(
                    identifier=str(last),
                    from_node_id=str(number),
                    to_node_id=str(marker["first"]),
                )
            elif coefficient > 0:
                reaction.add_segment(
                    identifier=str(last),
                    from_node_id=str(number),
                    to_node_id=str(marker["last"]),
                )
            else:
                # TODO: add exception
                raise Warning(f'Coefficient of "{metabolite}" cannot be 0')
        # Verify the number of Segments. They cannot have the same key of other
        # Segments from other reactions.

    def add_reaction(
        self,
        row: int,
        column: int,
        string: str,
        name: str,
        identifier: str,
        vertical: bool,
    ):
        """
        Parses and add given reaction string as a reaction for the
        JsonDictionary. It will automatically create all the necessary nodes
        and segments for the JSON.

        Args:
            string (str): Reaction string to be parsed.
            identifier (str): Identifier for the reaction.
            row (int): Row number from the visualization matrix.
            column (int): Column number of the visualization matrix.
            name (str): The name of the reaction.
            vertical (bool): If reaction should be displayed vertically.
        """
        # Add general data
        self._overview[identifier] = {
            "reactants": {},
            "products": {},
            "nodes": {},
            "position": Position(row=row, column=column),
            "index": "",
        }
        # Extract information for new reaction, nr of metabolites (string
        # representation)
        metabolite_dict = _convert_string(string=string)
        left_edge = self.R_WIDTH * column
        top_edge = self.R_HEIGHT * row
        # Labels
        if vertical:
            label_x = left_edge + (self.R_WIDTH / 2) + 40
            label_y = top_edge + (self.R_HEIGHT / 2)
        else:
            label_x = (left_edge + self.R_WIDTH / 2) - len(identifier) / 2 * 18
            label_y = top_edge + (self.R_HEIGHT) / 10
        # TODO: change this part
        reversibility = True
        reaction = Reaction(
            name=name,
            bigg_id=identifier,
            reversibility=reversibility,
            label_x=label_x,
            label_y=label_y,
            gene_reaction_rule="",
            genes=[],
            segments=dict(),
        )
        # Add nodes (metabolites and markers)
        self.map_metabolites(
            metabolite_dict=metabolite_dict,
            reaction=reaction,
            top_edge=top_edge,
            left_edge=left_edge,
            vertical=vertical,
        )
        self._add_reaction_markers(
            identifier=identifier,
            left_edge=left_edge,
            top_edge=top_edge,
            vertical=vertical,
        )
        # Add visual segments to reaction
        self.add_segments(reaction=reaction, metabolite_dict=metabolite_dict)
        # Define reaction number
        number = self._get_last_number(item="reactions")
        self.data["reactions"].update({str(number): reaction})
        self._overview[identifier]["index"] = str(number)
        debug_log.info(f'Reaction "{identifier}" added to the JsonDictionary.')

    def _reset(self):
        """
        Resets the data of JsonDictionary to its initial state. This helps
        when calling the function multiple times.
        """
        self._overview = {}
        self.data["reactions"] = {}
        self.data["nodes"] = {}
        self.reaction_scale = dict()

    def color_grading(
        self,
        color: List[Union[str, List[int], None]],
        min_max: List[float] = None,
        quantile: bool = False,
        max_steps: int = 100,
        n_steps: int = None,
    ):

        """
        Function that creates a color scale between two predefined colors. The
        number of color gradations corresponds to the number of fluxes.

        Parameters:
            color (list of str or list of list of int):
                list of two colors. The first color defines the endpoint for
                positive values, the second for negative values. The colors
                must be passed in their rgb representation.
            min_max(list of float,optional):
                List consisting of two floats. The first int describes the
                minimum, the second the maximum. These two values are
                additionally added as data values and at the same time ensure
                that all values greater or less than these are ignored when
                creating the gradient.
            quantile(bool, optional):
                Defines whether quantiles are to be used for the creation of
                the gradient. Otherwise, the steps are equally distributed
                between the minimum and maximum.
            n_steps (int, optional):
                Sets the number of color steps.
            max_steps(int, optional):
                Sets the maximum number of color steps.
        """

        # check if any flux values exist otherwise return
        if self.flux_solution is None:
            return

        # turn int arrays into numpy arrays
        color_positive = _color2np_rgb(color[0])
        color_negative = _color2np_rgb(color[1])
        color_intermediate = np.array([220, 220, 220], dtype=np.float32)

        try:
            flux = list(self.flux_solution.values())
        except AttributeError:
            flux = []

        # divide positive and negative values
        # bools used to handle the situation when both values of min_max are
        # positive or negative
        positive, negative, both_positive, both_negative = _divide_values(
            flux=flux, min_max=min_max
        )

        # array that will contain the configuration for escher
        reaction_scale = []

        # necessary variables for the loop that processes the positive
        # variables
        steps = len(positive)
        steps = min(steps, max_steps)

        # manually set steps overwrite the calculated steps
        if n_steps is not None:
            steps = math.floor(n_steps / 2)

        flux_values: ndarray

        try:
            if quantile:
                flux_values = np.linspace(0.0, 1.0, steps)

                # throws IndexError if positive == []
                flux_values = np.quantile(positive, flux_values)
            else:
                maximum = max(positive)
                start_value = maximum / steps

                # if both min_max values are positive the start is shifted
                # from zero to the set minimum
                if both_positive:
                    assert min_max is not None
                    start_value = min_max[0]

                flux_values = np.linspace(start_value, maximum, steps)

            with np.errstate(divide="ignore"):
                step_color = (color_intermediate - color_positive) / steps
            flux_color = color_intermediate.copy() - step_color

            flux_values.sort()
        except (ValueError, IndexError):
            flux_values = np.empty(shape=(0, 0))

        for flux_value in flux_values:
            reaction_scale.append(
                {
                    "type": "value",
                    "value": flux_value,
                    "color": "rgb(%d,%d,%d)"
                    % (flux_color[0], flux_color[1], flux_color[2]),
                }
            )

            flux_color -= step_color

        # add the intermediate step
        if not both_positive and not both_negative:
            reaction_scale.append(
                {"type": "value", "value": 0, "color": "rgb(220,220,220)"}
            )

        # updating of the variables for the negative portion of the flux values
        steps = len(negative)
        steps = min(steps, max_steps)

        # manually set steps overwrite the calculated steps
        if n_steps is not None:
            steps = math.floor(n_steps / 2)

        try:
            if quantile:
                flux_values = np.linspace(0.0, 1.0, steps)

                # throws IndexError if negative == []
                flux_values = np.quantile(negative, flux_values)
            else:
                minimum = min(negative)
                start_value = minimum / steps

                # if both min_max values are negative the start will be
                # shifted from zero to the set maximum
                if both_negative:
                    assert min_max is not None
                    start_value = min_max[1]
                flux_values = np.linspace(start_value, minimum, steps)

            with np.errstate(divide="ignore"):
                step_color = (color_intermediate - color_negative) / steps
            flux_color = color_intermediate - step_color

            flux_values.sort()
        except (ValueError, IndexError):
            flux_values = np.empty(shape=(0, 0))

        for flux_value in flux_values[::-1]:
            reaction_scale.append(
                {
                    "type": "value",
                    "value": flux_value,
                    "color": "rgb(%d,%d,%d)"
                    % (flux_color[0], flux_color[1], flux_color[2]),
                }
            )

            flux_color -= step_color

        self.reaction_scale = reaction_scale

    def visualize(
        self,
        filepath: Optional[Path] = None,
        vertical: bool = False,
        color: Optional[List[Union[str, List[int], None]]] = None,
        min_max: Optional[List[float]] = None,
        quantile: bool = False,
        max_steps: int = 100,
        n_steps: Optional[int] = None,
    ):
        """
        Saves the visualization of the JsonDictionary in given path as a HTML.
        Returns the builder for the JsonDictionary. If method is called in
        Jupyter or Qtconsole, it will show the embedded builder of the escher
        visualization. Else, it will open the default browser of the operating
        system and will load the previously saved HTML.

        .. note::
        Blank spaces are removed from the reactions.

        Args:
            filepath : Path, optional
                Path for the HTML. Defaults to "pathway.html" in the current
                working directory
            vertical : bool, optional
                Defines if pathway should be illustrated vertically.
                Defaults to False.
            color : list of str or list of lost of int, optional
                A list consisting of two entries. These define the endpoints
                for a color gradient. The entries can either be the colors as
                a list with three elements that define the RGB values or
                a string that defines the color name according to the css
                standard. If None is used here, no gradient will be generated.
                Defaults to None.
                Example: ["orange", "green"] or [[255, 165, 0], [0, 128, 0]]
            min_max: list of float, optional
                List consisting of two ints. The first int describes the
                minimum, the second the maximum. These two values are
                additionally added as data values and at the same time ensure
                that all values greater or less than these are ignored when
                creating the gradient.
                With Quantile = False, a data-independent gradient is created.
            quantile : bool, optional
                Defines whether quantiles are to be used for the creation of
                the gradient. Otherwise, the steps are equally distributed
                between the minimum and maximum.
            max_steps : int, optional
                Sets the maximum number of color steps.
            n_steps : int, optional
                Sets the maximum number of color steps.

        Returns:
            Builder: Escher builder object for the visualization

        See Also:
            Color names according to the css standard:
            https://www.w3schools.com/cssref/css_colors.asp
        """
        # Define path
        if not filepath:
            filepath = Path.cwd().joinpath("pathway.html")
        # Use relationship
        mapping = get_mapping(graph=self.graph)
        if vertical:
            mapping = transpose(matrix=mapping)
            # Modify Reaction-Box
            self.R_HEIGHT, self.R_WIDTH = self.R_WIDTH, self.R_HEIGHT
        # Modify canvas
        self.CANVAS_HEIGHT = self.R_HEIGHT * len(mapping)
        self.CANVAS_WIDTH = self.R_WIDTH * len(mapping[0])
        # Use reaction information
        for index_j, row in enumerate(mapping):
            for index_i, reaction in enumerate(row):
                # Add reactions only not 0
                if reaction == 0:
                    continue
                self.add_reaction(
                    row=index_j,
                    column=index_i,
                    name=reaction,
                    string=self.reaction_strings[reaction],
                    identifier=reaction,
                    vertical=vertical,
                )

        if color is not None:
            self.color_grading(
                color=color,
                min_max=min_max,
                quantile=quantile,
                max_steps=max_steps,
                n_steps=n_steps,
            )

        builder = Builder(
            # Check how reaction_styles behaves
            reaction_styles=["color", "text"],
            map_name=self.data["head"]["map_name"],
            map_json=self.json_dump(),
            reaction_scale=self.reaction_scale,
        )
        # This statement is needed, otherwise, all reactions labels will
        # appear with "(nd)".
        if self.flux_solution:
            builder.reaction_data = self.flux_solution
        builder.save_html(filepath=filepath)
        # builder.reaction_styles = ["color"]
        debug_log.info(f'Visualization located in "{filepath}"')
        # If in Jupyter, launch embedded widget. Otherwise, launch web-browser
        if not _in_notebook():
            # The context manager removes the ResourceWarning
            with catch_warnings():
                simplefilter(action="ignore", category=ResourceWarning)
                web_open("file://" + str(filepath))
        # Cleaning Up
        self._reset()
        return builder
