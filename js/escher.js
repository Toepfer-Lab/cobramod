    import {Builder} from "escher";

    function render({ model, el }) {

        let elem = document.createElement("div");
        el.appendChild(elem);

        let reaction_styles = model.get("reaction_styles");
        let map_name = model.get("map_name");
        let map_json = model.get("map_json");
        let reaction_data = model.get("reaction_data");
        let reaction_scale = model.get("reaction_scale");

        window.builder = Builder(
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