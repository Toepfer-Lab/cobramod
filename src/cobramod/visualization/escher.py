"""
This module contains an alternative Python integration for `Escher <https://github.com/zakandrewking/escher>`_ .
"""
from importlib import resources
from typing import Any, Optional, Dict, TypedDict, Literal, List

import anywidget
from traitlets import traitlets

from cobramod import static


class ReactionScale(TypedDict):
    """
    A Python class that represents the options available in Escher. It defines one point on the color scale.

    Args:
        type: The type of definition. A specific value can be used here (value), by means of
            quantiles (“Q1” or “Q3”) or by means of minimum maximum median or mean.
        value: It is only used to set the value if type value was used.
        color: The color to be used for the specified area.
        size: The width of the reactions covered by this definition.

    Notes:
        ReactionScale should only be an orientation for the possible options in Python. For further explanations and
        examples, please refer to `Escher's JavaScript documentation for 'escher.Builder.options.reaction_scale'. <https://escher.readthedocs.io/en/latest/javascript_api.html#escher.Builder.options.reaction_scale>`_
    """
    type: Literal["min", "max", "mean", "median", "Q1", "Q3", "value"]
    value: Optional[float]
    color: str
    size: int


class EscherIntegration(anywidget.AnyWidget):
    """
    An alternative Python integration for `Escher <https://github.com/zakandrewking/escher>`_ .

    """

    def __init__(
            self,
            map_name: Optional[str] = None,
            map_json: Optional[str] = None,
            reaction_data: Optional[Dict[str, float]] = None,
            reaction_scale: Optional[List[ReactionScale]] = None,
            reaction_styles=None,
            never_ask_before_quit: bool = False,
            *args: Any,
            **kwargs: Any,
    ):
        """

        Args:
            map_name: The name to be used by Escher for the map.
            map_json: The data structure to be displayed is encoded in JSON. Please refer to
                `Escher's documentation <https://escher.readthedocs.io/en/latest/convert_maps.html>`_
                for more information or to the
                `JSON schema <https://github.com/zakandrewking/escher/blob/master/jsonschema/1-0-0>`_.
            reaction_data: A dictionary of key-value pairs, where the keys are reaction IDs and the values are flux
                values.
            reaction_scale: A list consisting of ReactionScale.
                This list must consist of at least two ReactionScales when it is used.
            reaction_styles:
            never_ask_before_quit:
            *args:
            **kwargs:
        """
        super().__init__(*args, **kwargs)

        self.reaction_styles = reaction_styles
        self.map_name = map_name
        self.map_json = map_json
        self.reaction_scale = reaction_scale
        self.reaction_data = reaction_data
        self.never_ask_before_quit = never_ask_before_quit

    reaction_styles = traitlets.List(allow_none=True).tag(sync=True)
    map_name = traitlets.Unicode(allow_none=True).tag(sync=True)
    map_json = traitlets.Unicode(allow_none=True).tag(sync=True)
    reaction_scale = traitlets.Any(allow_none=True).tag(sync=True)
    reaction_data = traitlets.Dict(allow_none=True).tag(sync=True)
    never_ask_before_quit = traitlets.Bool(allow_none=False).tag(sync=True)

    _esm = resources.read_text(static, "escher.mjs")
