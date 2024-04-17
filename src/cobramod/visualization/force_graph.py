"""
This module contains the logic to create a three-dimensional representation
from a :py:class:`cobra.core.Group` or a :py:class:`cobra.Reaction`
using `3d-force-graph <https://github.com/vasturiano/3d-force-graph>`_ .
"""

import tempfile
import webbrowser
from dataclasses import dataclass, field
from typing import Union, Literal, Type, Optional

from cobra import Metabolite, Reaction, Solution
from cobra.core import Group


@dataclass(frozen= True)
class Nodes:
    id: str
    group: Literal["metabolite", "reaction"]

    def to_json(self) -> str:
        return f"""{{"id":"{self.id}", "group":"{self.group}"}}"""


@dataclass(unsafe_hash = True)
class Links:
    source: str = field(hash = True, compare = True)
    target: str = field(hash = True, compare = True)
    value: float = field(hash = False, compare = True)

    def to_json(self) -> str:
        return f"""{{"source":"{self.source}","target":"{self.target}","value":{self.value}}}"""


@dataclass()
class GraphData:
    nodes: set[Optional[Nodes]] = field(default_factory= set)
    links: set[Optional[Links]] = field(default_factory= set)

    def to_json(self) -> str:

        node_str = ",".join(x.to_json() for x in self.nodes)
        link_str = ",".join(x.to_json() for x in self.links)

        return f"""{{"nodes": [{node_str}],"links": [{link_str}]}}"""


def _reac2dict(reaction: Reaction, flux: float = 1) -> GraphData:
    nodes: set[Nodes] = set()
    links: set[Links] = set()

    nodes.add(
        Nodes(
            id = reaction.id,
            group= "reaction"
        )
    )

    for metabolite in reaction.metabolites:
        nodes.add(
            Nodes(
                id = metabolite.id,
                group = "metabolite"
            )
        )

        stoichiometry = reaction.get_coefficient(metabolite_id= metabolite.id)

        if stoichiometry * flux < 0:
            link:Links = Links(
                source= metabolite.id,
                target= reaction.id,
                value= - stoichiometry * flux
            )

        else:
            link:Links = Links(
                source= reaction.id,
                target= metabolite.id,
                value= stoichiometry * flux
            )

        links.add(link)

    data:GraphData = GraphData(
        nodes = nodes,
        links= links
    )

    return data

def _group2dict(group:Group, solution: Solution = None) -> GraphData:

    data : GraphData = GraphData()

    if solution is not None:
        if isinstance(solution, Solution):
            solution = solution.fluxes.to_dict()

    for member in group.members:

        if type(member) is Group:
            result: GraphData = _group2dict(group= member, solution= solution)

        elif type(member) is Reaction:
            if solution is not None:
                flux = solution.get(member.id, 1)
            else:
                flux = 1

            result: GraphData = _reac2dict(reaction= member, flux = flux)

        elif type(member) is Metabolite:
            node = set()
            node.add(Nodes(id= member.id, group="metabolite"))

            result:GraphData = GraphData(nodes= node)
        else:
            raise TypeError

        data.nodes.update(result.nodes)
        data.links.update(result.links)

    return data


class ForceGraphBuilder:
    def __init__(self, obj: Union[Type[Group], Type[Reaction]], solution: Solution = None):
        self.data: GraphData

        if isinstance(obj, Group):
            self.data = _group2dict(obj, solution = solution)

        elif isinstance(obj, Reaction):
            self.data = _reac2dict(obj)

        else:
            raise TypeError

    def _html(self) -> str:

        html =f"""
        <!DOCTYPE html>
        <html>
            <head>
                <script src="https://unpkg.com/3d-force-graph@1.73.3/dist/3d-force-graph.min.js"></script>
            </head>
            <body>
            <div id="3d-graph">
            </div>
            <script>
                let elem = document.getElementById("3d-graph")
                let graph_data = {self.data.to_json()}
                
                var Graph = ForceGraph3D()(elem)
                  .graphData(graph_data)
                  .nodeLabel("id")
                  .linkOpacity(1)
                  .linkAutoColorBy("value")
                  .linkDirectionalParticles(1)
                  .linkDirectionalParticleSpeed(d => d.value * 0.001)
                  .linkDirectionalParticleWidth(4)
                  .warmupTicks(100)
                  .cooldownTicks(0)
    
            </script>
            </body>
        </html>
        """

        return html

    def _repr_html_(self):
        return """
        <script>
        require.config({paths: {"3d-force-graph": "https://unpkg.com/3d-force-graph@1.73.3/dist/3d-force-graph.min.js"}});
        console.log(3d-force-graph)
        require(["3d-force-graph"], function(3d-force-graph) {
          console.log(3d-force-graph.version);
        });
        </script>
        """

    def open_html(self):
        html = self._html()

        tmp = tempfile.NamedTemporaryFile(suffix=".html")
        with open(tmp.name, 'w') as f:
            f.write(html)

        from IPython.display import IFrame
        IFrame(tmp.name, width=500, height=1000)
        webbrowser.open(tmp.name, new=0, autoraise=True)

