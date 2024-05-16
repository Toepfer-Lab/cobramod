"""
.. versionadded:: 1.3.0

This module contains an alternative Python integration for `Escher <https://github.com/zakandrewking/escher>`_ .
"""

from importlib import resources
from pathlib import Path
from typing import Any, Optional, Dict, TypedDict, Literal, List, Union

import anywidget
from traitlets import traitlets

from cobramod import static


class ReactionScale(TypedDict):
    """
    .. versionadded:: 1.3.0

    A Python class that represents the options available in Escher. It defines one point on the color scale.

    Args:
        type (Literal["min", "max", "mean", "median", "Q1", "Q3", "value"]): The type of definition. A specific value can be used here (value), by means of
            quantiles (“Q1” or “Q3”) or by means of minimum maximum median or mean.
        value (Optional[float]): It is only used to set the value if type value was used.
        color (str): The color to be used for the specified area.
        size (int): The width of the reactions covered by this definition.

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
    .. versionadded:: 1.3.0

    An alternative Python integration for `Escher <https://github.com/zakandrewking/escher>`_ .

    """

    def __init__(
        self,
        map_name: Optional[str] = None,
        map_json: Optional[str] = None,
        reaction_data: Optional[Dict[str, float]] = None,
        reaction_scale: Optional[List[ReactionScale]] = None,
        reaction_styles: Optional[List[Any]] = None,
        never_ask_before_quit: bool = False,
    ):
        """

        Args:
            map_name: The name to be used by Escher for the map.
            map_json: The data structure to be displayed encoded in JSON. Please refer to
                `Escher's documentation <https://escher.readthedocs.io/en/latest/convert_maps.html>`_
                for more information or to the
                `JSON schema <https://github.com/zakandrewking/escher/blob/master/jsonschema/1-0-0>`_.
            reaction_data: A dictionary of key-value pairs, where the keys are reaction IDs and the values are flux
                values.
            reaction_scale: A list consisting of :py:class:`~cobramod.visualization.escher.ReactionScale`.
                This list must consist of at least two ReactionScales when it is used.
            reaction_styles: Style options for the reactions see `Escher's JavaScript API <https://escher.readthedocs.io/en/latest/javascript_api.html#escher.Builder.options.reaction_styles>`_ for options and formatting.
            never_ask_before_quit: Option to control whether a warning dialog is displayed when the Escher Builder window is closed. See `Escher's JavaScript API <https://escher.readthedocs.io/en/latest/javascript_api.html#escher.Builder.options.never_ask_before_quit>`_ for reference.
        """
        super().__init__()

        self.map_name = map_name
        self.map_json = map_json
        self.reaction_scale = reaction_scale
        self.reaction_data = reaction_data
        self.reaction_styles = reaction_styles
        self.never_ask_before_quit = never_ask_before_quit

    reaction_styles: Optional[List[Any]] = traitlets.List(allow_none=True).tag(
        sync=True
    )  # type: ignore
    map_name = traitlets.Unicode(allow_none=True).tag(sync=True)
    map_json = traitlets.Unicode(allow_none=True).tag(sync=True)
    reaction_scale = traitlets.Any(allow_none=True).tag(sync=True)
    reaction_data: Optional[Dict[str, float]] = traitlets.Dict(
        allow_none=True
    ).tag(sync=True)  # type: ignore
    never_ask_before_quit = traitlets.Bool(allow_none=False).tag(sync=True)

    _esm = resources.read_text(static, "escher.mjs")

    @property
    def model_data(self):
        return "null"

    @property
    def embedded_css(self):
        return "null"

    @property
    def options(self):
        return {
            "reaction_styles": self.reaction_styles,
            "reaction_data": self.reaction_data,
            "reaction_scale": self.reaction_scale,
            "never_ask_before_quit": "true"
            if self.never_ask_before_quit
            else "false",
        }

    def save_html(self, filepath: Union[str, Path]):
        """
        This method creates a standalone HTML file that contains all the data of the
        Escher map and loads it automatically when it is opened.

        Args:
            filepath: The file in which the HTML file is to be saved.

        """

        html_rep = f"""
        <!DOCTYPE html>
        <html lang="en">
          <head>
            <title>Escher</title>
        
            <script src="https://unpkg.com/escher@1.7.3/dist/escher.min.js"></script>
        
            <meta charset="utf-8"/>
            <meta name="viewport" content="width=device-width, height=device-height,initial-scale=1.0, maximum-scale=1.0, user-scalable=no, minimal-ui"/>
          </head>
          <body>
            <div style="height: 100%; width: 100%;" id="map-container"></div>
        
            <script>
             escher.Builder(
                {self.map_json}, 
                {self.model_data}, 
                {self.embedded_css},
                escher.libs.d3_select('#map-container'), 
                {self.options},
             );
            </script>
          </body>
        </html>
        """
        if isinstance(filepath, str):
            filepath = Path(filepath)

        with open(filepath, "w") as f:
            f.write(html_rep)
