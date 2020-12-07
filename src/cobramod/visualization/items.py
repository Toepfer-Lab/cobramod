from contextlib import suppress
from collections import UserDict
from typing import Any, Dict, List

from cobramod.error import NodeAttributeError


class Node(UserDict):
    """
    Simple class that represent a node for the JSON schema of Escher. The
    node type can be either a 'metabolite', 'midmarker' or a 'multimarker'.
    Markers only need the arguments x and y. Check method
    :func:`cobramod.visualization.items.Node.__init__` to see all the
    corresponding arguments.
    """

    def __init__(
        self: Any,
        node_type: str,
        x: float,
        y: float,
        label_x: float = None,
        label_y: float = None,
        bigg_id: str = "",
        name: str = "",
        node_is_primary: bool = False,
    ):
        """
        Build a node for JsonCobramod, which can be parsed later into a proper
        JSON for Escher. A node can be a metabolite, midmarker or multimarker.
        Markers are the visible dots that are located in the middle of the
        reaction that defines the orientation (vertical, horizonzal).
        Metabolites represent the orange dots.

        Args:
            node_type: Can either be a 'metabolite', 'midmarker' or
                'multimarker'. Metabolites must have all arguments. Markers
                only need y and x.
            x (float): Location in x-axis for the node.
            y (float): Location in y-axis for the node
            label_x (float): Location in x-axis for the label of the node.
            label_y (float): Location in x-axis for the label of the node.
            bigg_id (str, optional): Identifier for the metabolite. It does not
                have to be a bigg identifier.
            name (str, optional): The name for the metabolite.
            node_is_primary (bool, optional): If the metabolite is a primary
                compound.
        """
        # Using a copy of locals, without self and the kwargs
        kwargs = {
            key: value
            for key, value in locals().copy().items()
            if key not in ("kwargs", "__class__", "self")
        }
        for method in (self.metabolite, self.marker):
            with suppress(NodeAttributeError):
                method(**kwargs)
        if "node_type" not in self.data.keys():
            raise NodeAttributeError(
                "Wrong attributes. Check 'node_types' and fix the arguments."
            )

    def marker(self, node_type: str, x: float, y: float, **kwargs):
        """
        Sets the attribute of the dictionary if the dictionary belongs to a
        midmarker or multimarker. Arguments specified in __init__.
        """
        # The kwargs is just to get rid of arguments that do not belong here
        # Must be yield, otherwise, the for loop wont work
        if node_type in ("midmarker", "multimarker"):
            super().__init__(self, node_type=node_type, x=x, y=y)
        else:
            raise NodeAttributeError(
                "Given node type does not represent a marker"
            )

    def metabolite(
        self,
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
        metabolite. Arguments specified in __init__,
        """
        # Must be yield, otherwise, the for loop wont work
        if node_type == "metabolite":
            kwargs = {
                key: value
                for key, value in locals().copy().items()
                if key not in ("kwargs", "__class__", "self")
            }
            super().__init__(self, **kwargs)
        else:
            raise NodeAttributeError(
                "Given node type does not represent a metabolite"
            )


class Segment(UserDict):
    """
    Build a segment for th reaction object for JsonCobramod, which can be
    later parsed into a proper JSON for Escher. A segment needs a starting
    position and and end-position. Additionally, the dots for the
    modification of the arrows (b1, b2) can be included.
    """

    def __init__(
        self,
        from_node_id: str,
        to_node_id: str,
        b1: Dict[str, float] = None,
        b2: Dict[str, float] = None,
    ):
        """
        Build a segment for the reactionn object that will be included in
        JsonCobramod. Check, documentation of the class for arguments of this
        method.

        Args:
            from_node_id (str): Identifier for the starting node.
            to_node_id (str): identifier for the end node.
            b1 (dict, optional): dictionary with x and y for the blue dot that
                shapes the arrow.
            b2 (dict, optional): dictionary with x and y for the blue dot that
            shapes the arrow.
        """
        kwargs = {
            key: value
            for key, value in locals().copy().items()
            if key not in ("kwargs", "__class__", "self")
        }
        super().__init__(self, **kwargs)


class Reaction(UserDict):
    """
    Simple class that represent a reaction for the JSON schema for Escher. A
    Reaction have the Segment that give the position of the reaction in the
    canvas. Check the method
    :func:`cobramod.visualization.items.Reaction.__init__' for the arguments
    for the creation of this special dictionary.
    """

    def __init__(
        self,
        name: str,
        bigg_id: str,
        reversibility: bool,
        label_y: float,
        label_x: float,
        gene_reaction_rule: str = "",
        genes: List[Dict[str, str]] = [],
        metabolites: List[Dict[str, str]] = [],
        segments: Dict[str, Segment] = {},
    ):
        """
        Creates an object with the information for the representation of a
        reaction in Escher.

        Args:
            name (str): Name of the reaction.
            bigg_id (str): Identifier of the reactions. It does not have to be
                from the database BIGG.
            reversibility (bool): If reaction is reversible
            label_y (float): Position in y-axis for the label.
            label_x (float): Position in x-axis for the label.
            gene_reaction_rule (str, optional): Gene rules that specify
                involved genes.
            genes (list, optional): A list with the genes involved for that
                identifier. Each item has to be a dictionary with the keys
                'bigg_id' and 'name'
            metabolites (list, optional): A list with the metabolites involved.
            Each iteam has to be an dictionary with the keys 'bigg_id' and
            'coefficient'.
            segments (dict, optional): dictionary with
                :func:`cobramod.visualization.items.Segment`.

        """
        kwargs = {
            key: value
            for key, value in locals().copy().items()
            if key not in ("kwargs", "__class__", "self")
        }
        super().__init__(self, **kwargs)
