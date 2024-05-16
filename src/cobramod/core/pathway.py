"""COBRApy Group-child class extension

The new class :class:`cobramod.pathway.Pathway` is child derived from
:class:`cobra.core.group.Group`. It extends some functionalities such as:

- solution: Obtain the solution for the specific members.
- visualize: get a :class:`escher.Builder` for that specific Pathway.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Optional, Union, Literal

import cobra.core as cobra_core

from cobramod.visualization.escher import EscherIntegration
from cobramod.visualization.force_graph import ForceGraphIntegration

try:
    import escher

    _has_escher = True
except ImportError:
    _has_escher = False

import pandas as pd

from cobramod.debug import debug_log
from cobramod.error import GraphKeyError
from cobramod.visualization.converter import JsonDictionary


class Pathway(cobra_core.Group):
    """
    A Sub-class from the original COBRApy :class:`cobra.Group`, which inherits
    all attributes and adds the method solution, to get a Solution for the
    members of this Class.

    Attributes
        vertical (bool, optional):
            Variable that determines whether the display should be vertical or
            horizontal using Escher. Defaults to False.
        color_negative (str or list of int, optional) :
            The color to use as the endpoint for the negative fluxes during the
            creation of the color gradient. All colors of the CSS standard can
            be used here or their RGB representation.
        color_positive (str or list of int, optional) :
            The color to use as the endpoint for the positive fluxes during the
            creation of the color gradient. All colors of the CSS standard can
            be used here or their RGB representation.
        color_min_max (list of float, optional) :
            The maximum and minimum to be taken into account when creating the
            color gradient. This creates these two values artificially to allow
            the creation of a data-independent color gradient. Fluxes larger or
            smaller are ignored accordingly.
        color_quantile (bool, optional) :
            Attribute that defines whether the color gradient should be
            determined through quantiles or equally distributed between the
            maximum and the minimum. Defaults to False which means that the
            gradations are evenly distributed.
        color_n_steps (int, optional) :
            The number of steps used when creating the color gradient. Uses the
            number of fluxes by default. The default value is None.
        color_max_steps (int, optional) :
            The maximum number of steps to use when creating the color
            gradient. Default value is 100.
    See Also:
            Color names according to the css standard:
            https://www.w3schools.com/cssref/css_colors.asp
    """

    def __init__(
        self,
        id: str,
        name: str = "",
        members: Optional[set[cobra_core.Object]] = None,
        kind: Optional[str] = None,
    ):
        """
        Regular __init__ of the class. It maintains the same attributes as its
        parent :class:`cobra.core.group.Group`. Included is the order of the
        reactions.

        Args:
            id (str): Identifier of the pathway.
            name (str): Name for the pathway.
            members (set): Includes the members for the pathway. It should
                include only :class:`cobra.core.reaction.Reaction`. Defaults to
                None.
            kind (str): The kind of group, as specified for the Groups feature
                in the SBML level 3 package specification. Can be any of
                "classification", "partonomy", or "collection". The default is
                "collection". Please consult the SBML level 3 package
                specification to ensure you are using the proper value for
                kind. (Text extracted from parent.)
        """
        super().__init__(
            id=id,
            name=name,
            kind=kind,
        )

        # Loop has to be after __init__, otherwise, behavior of class changes.
        self.graph = dict()  # type: ignore
        self.notes: dict[str, Any] = {"ORDER": dict()}

        if members:
            for member in members:
                # Remove no Reactions
                if not isinstance(member, cobra_core.Reaction):
                    debug_log.debug(
                        f'Member "{member.id}" is not a Reaction. This object will'
                        f'not be included in the pathway "{self.id}".'
                    )
                    continue
                super().add_members(new_members=[member])

        # Attributes that can be used for the customization of the
        # visualization
        self.color_negative: Optional[Union[str, list[int]]] = None
        self.color_positive: Optional[Union[str, list[int]]] = None
        self.color_min_max: Optional[list[float]] = None
        self.color_max_steps: int = 100
        self.color_n_steps: Optional[int] = None
        self.color_quantile: bool = False
        self.vertical: bool = False

    def _filter(
        self, solution: cobra_core.Solution, attribute: str
    ) -> pd.Series:
        """
        Filters given attribute of Solution. Only members of the Pathway will
        appear.
        """
        return getattr(solution, attribute).filter(
            items=[member.id for member in self.members], axis="index"
        )

    def add_members(self, new_members: list[cobra_core.Reaction]):
        """
        Add given list of :class:`cobra.core.reaction.Reaction` into the
        Pathway.

        Args:
            new_members (list): List of Reactions to add to the class.

        Raises:
            TypeError: If not all members are proper Reaction objects.
        """
        if not all(
            [isinstance(member, cobra_core.Reaction) for member in new_members]
        ):
            raise TypeError("Not all given members are Reactions.")

        super().add_members(new_members=new_members)
        for member in new_members:
            try:
                # There is no need to extend the tuple or strings because
                # the father-reaction has at least 1 reaction
                self.notes["ORDER"][member.id]
            except KeyError:
                self.notes["ORDER"][member.id] = None

    def solution(self, solution: cobra_core.Solution) -> cobra_core.Solution:
        """
        Returns a :class:`cobra.Solution` with only the members of the pathway.

        Args:
            solution (Solution): Original COBRApy :class:`cobra.Solution` to
                filter.
        Returns:
            Solution: Filtered solution containing only members of the Pathway
                class.
        """
        return cobra_core.Solution(
            objective_value=solution.objective_value,
            status=solution.status,
            fluxes=self._filter(solution=solution, attribute="fluxes"),
            reduced_costs=self._filter(
                solution=solution, attribute="reduced_costs"
            ),
            shadow_prices=self._filter(
                solution=solution, attribute="shadow_prices"
            ),
        )

    @classmethod
    def _transform(cls, obj: cobra_core.Group) -> "Pathway":
        """
        Transform given :class:`cobra.Group` into a proper cobramod Pathway
        object.
        """
        pathway = cls(
            id=obj.id, name=obj.name, members=obj.members, kind=obj.kind
        )
        pathway._model = obj._model
        # Add directly from the notes or create a quick order
        try:
            # When loading model, COBRApy will modify the string
            if isinstance(obj.notes["ORDER"], str):
                obj.notes["ORDER"] = eval(
                    obj.notes["ORDER"].replace("&apos;", '"')
                )
            pathway.graph = pathway.notes["ORDER"] = obj.notes["ORDER"]

        except KeyError:
            # Check only Reactions
            reactions = [
                rxn
                for rxn in obj.members
                if isinstance(rxn, cobra_core.Reaction)
            ]
            order = {}
            for value, reaction in enumerate(reactions):
                try:
                    order[reaction.id] = reactions[value + 1].id
                except IndexError:
                    # It should be the last
                    order[reaction.id] = None
            pathway.graph = pathway.notes["ORDER"] = order

        debug_log.info(
            f'Group-object "{pathway.id}" was transformed to a Pathway-object.'
        )
        return pathway

    def modify_graph(self, reaction: str, next_reaction: Union[str, None]):
        """
        Modifies the order of the graph. This is useful when merging multiple
        pathways or joining reactions. In the graph, the selected reaction
        will be forced to show "next_reaction" as its successor.

        Args:
            reaction (str): Identifier of the reaction to modify in the graph.
            next_reaction (str, None): Identifier of the next reaction. This
                reaction will take place after "reaction". If None is passed,
                then "reaction" will not have successors.

        Raises:
            GraphKeyError: If the reaction or the next_reaction does not appear
            in the graph of the pathway.
        """
        # FIXME: add behavior for changing not NoneTypes
        # Pathway.graph is responsable for the mapping
        if next_reaction not in self.graph.keys() and next_reaction:
            raise GraphKeyError(
                f'Pathway "{self.id}" does not have a reaction '
                + f'{next_reaction}". Check that the reaction exist.'
            )
        try:
            self.graph[reaction]
            self.graph[reaction] = next_reaction
            # We need to modify the notes as well because it will change the
            # notes
            msg = f'The reaction order of pathway "{self.id}" was modified. '
            if not next_reaction:
                msg += f'There is no next reaction after "{reaction}".'
            else:
                msg += (
                    f'Reaction "{next_reaction}" takes place after '
                    + f'"{reaction}".'
                )
            self.notes["ORDER"][reaction] = next_reaction
            debug_log.info(msg=msg)
        except KeyError:
            raise GraphKeyError(
                f'Pathway "{self.id}" does not have a reaction '
                + f'{reaction}". Check that the reaction exist.'
            )

    def visualize(
        self,
        solution_fluxes: Optional[
            Union[cobra_core.Solution, dict[str, float]]
        ] = None,
        filename: Optional[Union[str, Path]] = None,
        vis: Literal["escher", "escher-custom", "3d-force"] = "escher",
        never_ask_before_quit: bool = False,
    ) -> Union[escher.Builder, EscherIntegration, ForceGraphIntegration, None]:
        """
        .. versionchanged:: 1.3.0
            The 'vis' parameter has been added. This allows one to choose between different visualization tools.

        Returns a :class:`escher.Builder`, which can be used to create visual
        representations of the pathway.

        :param solution_fluxes: Series or Dictionary with fluxes. The values will be then showed in the Builder. Defaults to None.

        :param filename: Path for the HTML. Defaults to "pathway.html" in the current working directory.

        :param vis:
            .. versionadded:: 1.3.0
            Parameter that determines the visualization tool used. It is possible to choose between the original
            Escher integration [escher], the one embedded in CobraMod [escher-custom] and a 3-dimensional
            force directed graph visualization [3d-force].

            .. deprecated:: 1.3.0
                The original python integration of Escher will be removed in a future version due to dependency
                conflicts with Jupyter. The integration embedded in CobraMod will take its place in the future.
                This can already be used by setting 'vis' to "escher-custom".

        :param never_ask_before_quit:
            .. versionadded:: 1.3.0

            Option to control whether a warning dialog is displayed when the Escher Builder window is closed.
            Only has an effect when using Escher for visualization.

        """

        if vis == "3d-force":
            widget = ForceGraphIntegration()
            widget.model = self
            widget.solution = solution_fluxes

            return widget

        json_dict = JsonDictionary()
        if filename is None:
            filename = "pathway.html"

        # Define solution. If None, nothing will be added. Either dict or
        # regular solution
        if solution_fluxes is not None:
            if isinstance(solution_fluxes, cobra_core.Solution):
                solution_fluxes = solution_fluxes.fluxes.to_dict()
            json_dict.flux_solution = solution_fluxes

        # Get graph and add to json_dict
        json_dict.graph = self.graph.copy()
        reactions: dict[str, str] = {m.id: m.reaction for m in self.members}
        json_dict.reaction_strings = reactions

        if vis == "escher-custom":
            builder = json_dict.visualize(
                filepath=filename,
                vertical=self.vertical,
                color=[self.color_positive, self.color_negative],
                min_max=self.color_min_max,
                quantile=self.color_quantile,
                max_steps=self.color_max_steps,
                n_steps=self.color_n_steps,
                custom_integration=True,
                never_ask_before_quit=never_ask_before_quit,
            )
            return builder

        if _has_escher:
            builder = json_dict.visualize(
                filepath=filename,
                vertical=self.vertical,
                color=[self.color_positive, self.color_negative],
                min_max=self.color_min_max,
                quantile=self.color_quantile,
                max_steps=self.color_max_steps,
                n_steps=self.color_n_steps,
            )
            debug_log.info(f'Visualization saved in "{filename}"')
            warnings.warn(
                "The use of Escher's own Python integration will be removed in a future CobraMod version to "
                "avoid installation problems due to dependency incompatibilities. Instead, you can use "
                "CobraMod's built-in Escher integration. To use it, simply pass 'escher-custom' as an "
                "argument to 'vis'.",
                category=DeprecationWarning,
                stacklevel=2,
            )

            return builder
        else:
            warnings.warn(
                "Package Escher was not found. No visualization in available. "
                "You can install this extra-dependency running "
                "'pip install cobramod[escher]'. Alternatively, you can use CobraMod's own integration of Escher. "
                "To use it, just specify 'escher-custom' for the 'vis' argument."
            )
            return None

    def _repr_html_(self):
        """
        Returns a HTML string with the attributes of the Pathway
        """
        return f"""
<table> <tbody> <tr> <td><strong>Pathway identifier</strong></td>
<td>{self.id}</td> </tr> <tr> <td><strong>Name</strong></td>
<td>{self.name}</td> </tr> <tr> <td><strong>Memory address</strong></td>
<td>0x0{id(self)}</td> </tr> <tr> <td><strong>Reactions involved</strong></td>
<td> <p>{", ".join([rxn.id for rxn in self.members])}</p> </td> </tr> <tr>
<td><strong>Genes involved<br /></strong></td> <td> <p>{", ".join([gene.id for
rxn in self.members for gene in rxn.genes])}</p> </td> </tr> <tr>
<td><strong>Visualization attributes</strong></td> <td> <ul> <li>vertical =
{self.vertical}</li> <li>color_negative = {self.color_negative}</li>
<li>color_positive = {self.color_positive}</li> <li>color_quantile =
{self.color_quantile}</li> </ul> </td> </tr> </tbody> </table> <p>&nbsp;</p>"""


def model_convert(model: cobra_core.Model):
    """
    Converts all Group objects in the given model to proper cobramod
    :class:`cobramod.pathway.Pathway`
    """
    for index, group in enumerate(iterable=model.groups):
        model.groups[index] = Pathway._transform(obj=group)
