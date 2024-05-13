    import {Builder} from "escher";

    function render({ model, el }: { model: DOMWidgetModel; el: HTMLElement; }) {

        let elem = document.createElement("div");
        let cell = el.getBoundingClientRect()

        el.appendChild(elem);

        let reaction_styles = model.get("reaction_styles");
        let map_name: string = model.get("map_name");
        let map_json: string = model.get("map_json");
        let reaction_data = model.get("reaction_data");
        let reaction_scale = model.get("reaction_scale");
        let never_ask_before_quit: boolean = model.get("never_ask_before_quit")

        elem.style.width = cell.width.toString()+"px"
        elem.style.height = cell.width.toString()+"px"

        let builder = Builder(
            JSON.parse(map_json),
            null,
            null,
            elem,
            {
                "reaction_styles": reaction_styles,
                "reaction_data": reaction_data,
                "reaction_scale": reaction_scale,
                "never_ask_before_quit": never_ask_before_quit,
            },
        )


        model.on("change:reaction_scale", () => {
             builder.options.reaction_scale = model.get("reaction_scale");
        });

    }

    export default { render };