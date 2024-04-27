"""
This module contains an alternative Python integration for `Escher <https://github.com/zakandrewking/escher>`_ .
"""
from importlib import resources
from pathlib import Path
from typing import Any

import anywidget
from traitlets import traitlets

from cobramod import static


class EscherIntegration(anywidget.AnyWidget):
    def __init__(
        self,
        reaction_styles=None,
        map_name=None,
        map_json=None,
        reaction_scale=None,
        reaction_data=None,
        never_ask_before_quit:bool = False,
        *args: Any,
        **kwargs: Any,
    ):
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
