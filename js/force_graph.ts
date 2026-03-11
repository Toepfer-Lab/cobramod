import ForceGraph3D from '3d-force-graph';

function render({model, el}: { model: DOMWidgetModel; el: HTMLElement; }) {

    let elem = document.createElement("div");
    el.appendChild(elem);
    let Graph: any;

    // wait until el is finished and has a size
    const ro = new ResizeObserver((entries) => {
        for (const entry of entries) {
            const width = entry.contentRect.width;
            if (width > 0) {
                if (!Graph) {
                    let graph_data = JSON.parse(model.get("_model_rep"))
                    Graph = new ForceGraph3D(elem)
                        .graphData(graph_data)
                        .nodeLabel("id")
                        .linkOpacity(1)
                        .linkAutoColorBy("value")
                        .linkDirectionalParticles(1)
                        .linkDirectionalParticleSpeed(d => d["value"] * 0.001)
                        .linkDirectionalParticleWidth(4)
                        .warmupTicks(100)
                        .cooldownTicks(0)
                        .width(width)
                        .height(width / 2)

                    model.on("change:_model_rep", () => {
                        Graph.graphData(JSON.parse(model.get("_model_rep")))
                    });

                    model.on("msg:custom", msg => {
                        switch (msg.type) {
                            case "create_layout":
                                let nodes = {}
                                Graph.graphData().nodes.forEach((n) => {
                                    nodes[n.id] = {"x": n.x, "y": n.y, "z": n.z};
                                });
                                model.send({type: "layout", positions: nodes});
                                break;
                            case "load_layout":
                                let positions = msg.positions;
                                Graph.graphData().nodes.forEach((n) => {
                                    let pos = positions[n.id];
                                    n.fx = pos.x;
                                    n.fy = pos.y;
                                    n.fz = pos.z;
                                });
                                Graph.cooldownTicks(1)
                                Graph.d3ReheatSimulation()
                                break;
                            default:
                                console.log(`Unknown ${msg}.`);
                        }
                    });
                } else {
                    Graph.width(width).height(width / 2);
                }
            }
        }
    });
    ro.observe(el);
}

export default {render};