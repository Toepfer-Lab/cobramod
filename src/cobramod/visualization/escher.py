"""
This module contains an alternative Python integration for `Escher <https://github.com/zakandrewking/escher>`_ .
"""
from typing import Any

import anywidget
from traitlets import traitlets


class EscherIntegration(anywidget.AnyWidget):
    def __init__(self,
                 reaction_styles = None,
                 map_name = None,
                 map_json = None,
                 reaction_scale = None,
                 reaction_data = None,
                 *args: Any,
                 **kwargs: Any
                 ):
        super().__init__(*args, **kwargs)

        self.reaction_styles = reaction_styles
        self.map_name = map_name
        self.map_json = map_json
        self.reaction_scale = reaction_scale
        self.reaction_data = reaction_data

    reaction_styles= traitlets.List(allow_none=True).tag(sync=True)
    map_name= traitlets.Unicode(allow_none=True).tag(sync=True)
    map_json= traitlets.Unicode(allow_none=True).tag(sync=True)
    reaction_scale= traitlets.Any(allow_none=True).tag(sync=True)
    reaction_data= traitlets.Dict(allow_none=True).tag(sync=True)

    _esm = """
    import * as escher_imp from "https://unpkg.com/escher@1.7.3/dist/escher.min";

    function render({ model, el }) {
    
        let elem = document.createElement("div");
        el.appendChild(elem);  
    
        let reaction_styles = model.get("reaction_styles");
        let map_name = model.get("map_name");
        let map_json = model.get("map_json");
        let reaction_data = model.get("reaction_data");  
        let reaction_scale = model.get("reaction_scale");
            
        window.builder = escher.Builder(
            JSON.parse(map_json),
            null,
            null,
            elem,
            {
                "reaction_styles": reaction_styles,
                "reaction_data": reaction_data,
                "reaction_scale": reaction_scale,
            },
        )
        
        
        model.on("change:reaction_scale", () => {
             window.builder.options.reaction_scale = model.get("reaction_scale");
        });
        
    }

    export default { render };
    """
