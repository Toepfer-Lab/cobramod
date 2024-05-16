"""
.. versionadded:: 1.3.0

This module contains the logic to create a three-dimensional representation
from a :py:class:`cobra.core.Group` or a :py:class:`cobra.Reaction`
using `3d-force-graph <https://github.com/vasturiano/3d-force-graph>`_ .
"""

import csv
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Union, Literal, Type, Optional, Any

import anywidget
from cobra import Metabolite, Reaction, Solution
from cobra.core import Group
from traitlets import traitlets

from cobramod import static


@dataclass(frozen=True)
class Nodes:
    id: str
    group: Literal["metabolite", "reaction"]

    def to_json(self) -> str:
        return f"""{{"id":"{self.id}", "group":"{self.group}"}}"""


@dataclass(unsafe_hash=True)
class Links:
    __source: str = field(hash=True, compare=True)
    __target: str = field(hash=True, compare=True)
    value: float = field(hash=False, compare=True)

    def __init__(self, source: str, target: str, value: float):
        self.__source = source
        self.__target = target
        self.value = value

    @property
    def source(self) -> str:
        return self.__source

    @source.setter
    def source(self, value):
        raise AttributeError("The source attribute is read only. ")

    @property
    def target(self) -> str:
        return self.__target

    @target.setter
    def target(self, value):
        raise AttributeError("The target attribute is read only. ")

    def to_json(self) -> str:
        return f"""{{"source":"{self.source}","target":"{self.target}","value":{self.value}}}"""


@dataclass()
class GraphData:
    nodes: set[Nodes] = field(default_factory=set)
    links: set[Links] = field(default_factory=set)

    def to_json(self) -> str:
        node_str: str = ""
        link_str: str = ""

        if len(self.nodes) > 0:
            node_str = ",".join(x.to_json() for x in self.nodes)

        if len(self.links) > 0:
            link_str = ",".join(x.to_json() for x in self.links)

        return f"""{{"nodes": [{node_str}],"links": [{link_str}]}}"""

    def to_dict(self) -> dict:
        return {"nodes": [self.nodes], "links": [self.links]}


def _reac2dict(reaction: Reaction, flux: float = 1) -> GraphData:
    """
    Function that generates a :py:class:`GraphData` object from a Cobra.Reaction.

    Args:
        reaction: A :py:class:`cobra.Reaction` that is to be converted into a :py:class:`GraphData` object.
        flux: The flux of this reaction.

    Returns: A :py:class:`GraphData` object that contains the definition for the passed reaction.

    """

    nodes: set[Nodes] = set()
    links: set[Links] = set()

    nodes.add(Nodes(id=reaction.id, group="reaction"))

    for metabolite in reaction.metabolites:
        nodes.add(Nodes(id=metabolite.id, group="metabolite"))

        stoichiometry = reaction.get_coefficient(metabolite_id=metabolite.id)
        link: Links

        if stoichiometry * flux < 0:
            link = Links(
                source=metabolite.id,
                target=reaction.id,
                value=-stoichiometry * flux,
            )

        else:
            link = Links(
                source=reaction.id,
                target=metabolite.id,
                value=stoichiometry * flux,
            )

        links.add(link)

    data: GraphData = GraphData(nodes=nodes, links=links)

    return data


def _group2dict(
    group: Group, solution: Union[Solution, dict] = None
) -> GraphData:
    """
    Function that generates a :py:class:`GraphData` object from a :py:class:`cobra.Reaction`.

    Args:
        group: A :py:class:`cobra.Reaction` that is to be converted into a :py:class:`GraphData` object.
        solution: The flux of this reaction.

    Returns:
        A :py:class:`GraphData` object that contains the definition for the passed reaction.

    """

    data: GraphData = GraphData()

    if solution is not None:
        if isinstance(solution, Solution):
            solution = solution.fluxes.to_dict()

    result: GraphData

    for member in group.members:
        if type(member) is Group:
            result = _group2dict(group=member, solution=solution)

        elif type(member) is Reaction:
            if solution is not None:
                flux = solution.get(member.id, 1)
            else:
                flux = 1

            result = _reac2dict(reaction=member, flux=flux)

        elif type(member) is Metabolite:
            node = set()
            node.add(Nodes(id=member.id, group="metabolite"))

            result = GraphData(nodes=node)
        else:
            raise TypeError

        data.nodes.update(result.nodes)
        data.links.update(result.links)

    return data


class ForceGraphIntegration(anywidget.AnyWidget):
    """
    .. versionadded:: 1.3.0

    Widget for displaying a :py:class:`cobra.core.group.Group` or :py:class:`cobra.Reaction` as a force directed graph.
    """

    def __init__(self):
        """ """
        super().__init__()
        self.on_msg(self._handle_custom_msg)

    _model: Union[Type[Group], Type[Reaction], None] = None
    _solution: Optional[dict] = None

    # _model_rep is the data basis for the widget
    # Changes update the front end
    _model_rep: str = """{"nodes": [], "links": []}"""
    _model_rep = traitlets.Unicode().tag(sync=True)  # type: ignore

    @property
    def model(self) -> Optional[Union[Type[Group], Type[Reaction]]]:
        """
        The Model to be represented. It can ether be a :py:class:`cobra.core.group.Group` or :py:class:`cobra.Reaction`.
        It is set to None upon initialization.
        """

        return self._model

    @model.setter
    def model(self, value: Union[Type[Group], Type[Reaction]]):
        self._model = value
        self._create_model_rep()

    @property
    def solution(self) -> Optional[Union[Type[Solution], dict[str, float]]]:
        """
        The flux values to be taken into account when creating the graph. These scale the stoichiometry of the
        respective reaction. If a reaction is not found in this object, a flux of 1 is assumed. The flux values can
        either be passed as cobra.Solution or as dict, where the key is the reaction ID and the value is the flux
        value. It is set to None upon initialization.

        """
        return self._solution

    @solution.setter
    def solution(self, value: Union[dict[str, float], Solution]):
        if isinstance(value, Solution):
            value = value.fluxes.to_dict()

        self._solution = value
        self._create_model_rep()

    def _create_model_rep(self):
        """
        Function to create a representation whenever a model and/or a solution are assigned.
        This representation is used as a data basis for 3d-force-graph.
        The representation is assigned to 'self._model_rep'.
        """

        if self._model is None:
            return

        if isinstance(self._model, Group):
            data = _group2dict(self._model, solution=self._solution)

        elif isinstance(self._model, Reaction):
            if self._solution is not None:
                flux = self._solution.get(self._model.id, 1)
            else:
                flux = 1

            data = _reac2dict(self._model, flux=flux)

        else:
            raise TypeError

        self._model_rep = data.to_json()

    def save_layout(self, file: Union[str, Path]):
        """
        Method to save the layout of the current display. The resulting file is a CSV containing the ID
        of the metabolite and its X, Y and Z position in the other columns.

        Args:
            file: The path where the layout is to be saved.

        """

        if isinstance(file, str):
            file = Path(file)

        self.__csv_path = file
        self.send({"type": "create_layout"})

    def load_layout(self, file: Union[str, Path]):
        """
        Method to restore the layout that was saved using
        :py:meth:`~cobramod.visualization.force_graph.ForceGraphIntegration.save_layout`.

        Args:
            file: Path to the file containing the saved layout.

        """

        if isinstance(file, str):
            file = Path(file)

        with open(file) as open_file:
            reader = csv.DictReader(open_file)
            positions = {}
            for row in reader:
                key = row["ID"]
                x = row["x"]
                y = row["y"]
                z = row["z"]

                positions[key] = {"x": x, "y": y, "z": z}

        self.send({"type": "load_layout", "positions": positions})

    def _handle_custom_msg(self, data: dict, buffers: Any):
        if data["type"] == "layout":
            positions: dict[str, dict[Literal["x", "y", "z"], float]] = data[
                "positions"
            ]
            fieldnames = ["ID", "x", "y", "z"]
            with open(self.__csv_path, "w") as file:
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                for key, value in positions.items():
                    writer.writerow(
                        {
                            "ID": key,
                            "x": value.get("x", 0),
                            "y": value.get("y", 0),
                            "z": value.get("z", 0),
                        }
                    )

    _esm = resources.files(static).joinpath("force_graph.mjs").read_text()
