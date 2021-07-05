#!/usr/bin/env python3
"""COBRApy Group-child class extension

The new class :class:`cobramod.pathway.Pathway" is child derived from
:class:`cobra.core.group.Group`. It extends some functionalities such as:

- solution: Obtain the solution for the specific members.
- visualize: get a :class:`escher.Builder` for that specific Pathway.
"""
from pathlib import Path
from typing import Any, Dict, List, Union, Optional

from escher import Builder
from cobra.core.group import Group
from cobra.core.model import Model
from cobra.core.solution import Solution
from cobra.core.reaction import Reaction
from cobra.core.dictlist import DictList
from pandas import Series

from cobramod.debug import debug_log
from cobramod.error import GraphKeyError
from cobramod.visualization.converter import JsonDictionary


class Pathway(Group):
    """
    A Sub-class from the original COBRApy :class:`cobra.Group`, which inherits
    all its attributes and adds the method solution, to get a Solution for the
    members of this Class.

    Attributes
        vertical (bool, optional):
            Variable that determines whether the display should be vertical or
            horizontal using escher. Defaults to False.
        color_negative (str or list of int, optional) :
            The color to use as the endpoint for the negative fluxes during the
            creation of the color gradient. All colors of the css standard can
            be used here or their rgb representation.
        color_positive (str or list of int, optional) :
            The color to use as the endpoint for the positive fluxes during the
            creation of the color gradient. All colors of the css standard can
            be used here or their rgb representation.
        color_min_max (list of float, optional) :
            The maximum and minimum to be taken into account when creating the
            color gradient. This creates these two values artificially to allow
            the creation of a data-independent color gradient. Fluxes larger or
            smaller are ignored accordingly.
        color_quantile (bool, optional) :
            Attributes that determines whether the color gradient should be
            determined through quantiles or equally distributed between the
            maximum and the minimum. Defaults to None which means that the
            gradations are evenly distributed. Default value is False.
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
        members: DictList = DictList(),
        kind: str = None,
    ):
        """
        Regular __init__ of the class. It maintains the same attributes as its
        parent :class:`cobra.core.group.Group`. Included is the order of the
        reactions.

        Args:
            id (str): Identifier of the pathway.
            name (str): Name for the patway.
            members (DictList): Includes the members for the pathway. It should
                include only :class:`cobra.core.reaction.Reaction`. Defaults to
                None.
            kind (str): The kind of group, as specified for the Groups feature
                in the SBML level 3 package specification. Can be any of
                "classification", "partonomy", or "collection". The default is
                "collection". Please consult the SBML level 3 package
                specification to ensure you are using the proper value for
                kind. (Text extracted from parent.)
        """
        super().__init__(id=id, name=name, kind=kind)
        # Loop has to be after __init__, otherwise, behavior of class changes.
        # TODO: Is order necessary? Maybe use notes["order"] directly
        self.order: List[str] = list()
        self.graph: dict = dict()
        self.notes: Dict[str, Any] = {"ORDER": dict()}
        for member in members:
            # Remove no Reactions
            if not isinstance(member, Reaction):
                debug_log.debug(
                    f'Member "{member.id}" is not a Reaction. Skipped'
                )
                continue
            self.add_members(new_members=[member])

        # Attributes that can be used for the customization of the
        # visualization
        self.color_negative: Optional[Union[str, List[int]]] = None
        self.color_positive: Optional[Union[str, List[int]]] = None
        self.color_min_max: Optional[List[float]] = None
        self.color_max_steps: int = 100
        self.color_n_steps: Optional[int] = None
        self.color_quantile: bool = False
        self.vertical: bool = False

    def _filter(self, solution: Solution, attribute: str) -> Series:
        """
        Filters given attribute of Solution. Only members of the Pathway will
        appear.
        """
        return getattr(solution, "fluxes").filter(
            items=[member.id for member in self.members], axis="index"
        )

    def __check_copy(self):
        """
        This method checks if the length of the members is zero and attribute
        "order" is larger that zero. This check is made is to check if method
        :func:`cobra.core.model.Model.copy` is called. For copies, order must
        reset.
        """
        if len(self.members) == 0 and len(self.order) > 0:
            self.order: List[str] = list()
            self.graph: dict = dict()
            debug_log.debug(
                f'Attribute order from pathway "{self.id}" reset. '
                f"Check if a method copy() from cobra.Model was called."
            )

    def add_members(self, new_members: List[Reaction]):
        """
        Add given list of :class:`cobra.core.reaction.Reaction` into the
        Pathway.

        Args:
            new_members (list): List of Reactions to add to the class.

        Raises:
            TypeError: If not all members are proper Reaction objects.
        """
        if not all([isinstance(member, Reaction) for member in new_members]):
            raise TypeError("Not all given members are Reactions.")

        self.__check_copy()
        super().add_members(new_members=new_members)
        for member in new_members:
            self.notes["ORDER"][member.id] = None
        # Extend order in order to use it later for the visualization.
        self.order.extend((member.id for member in new_members))

    def solution(self, solution: Solution) -> Solution:
        """
        Returns a :class:`cobra.Solution` with only the members of the pathway.

        Args:
            solution (Solution): Original COBRApy :class:`cobra.Solution` to
                filter.
        Returns:
            Solution: Filtered solution with only members of the Pathway class.
        """
        return Solution(
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
    def _transform(cls, obj: Group) -> "Pathway":
        """
        Transform given :class:`cobra.Group` into a proper cobramod Pathway
        object.
        """
        return cls(
            id=obj.id, name=obj.name, members=obj.members, kind=obj.kind
        )

    def modify_graph(self, reaction: str, next_reaction: Union[str, None]):
        """
        Modifies the order of the graph. This is useful when merging multiple
        pathways or joining reactions. In the graph, the selected reaction
        will be forced to show "next_reaction" as  its successor.

        Args:
            reaction (str): Identifier of the reaction to modify in the graph
            next_reaction (str, None): Identifier of the next reaction. This
                reaction will take place after "reaction". If None is passed,
                then "reaction" will not have successors.

        Raises:
            GraphKeyError: If the reaction or the next_reaction does not appear
            in the graph of the pathway
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
        solution_fluxes: Union[Solution, Dict[str, float]] = None,
        filename: Path = None,
    ) -> Builder:
        """
        Returns a :class:`escher.Builder`, which can be use to create visual
        representations of the pathway.

        Args:
            solution_fluxes (Solution, dict): Series or Dictionary with fluxes.
                The values will be then showed in the Builder.
                Defaults to None.
            filename (Path): Path for the HTML. If None is passed, then
                default to "pathway.html" in the current working directory.
        """
        json_dict = JsonDictionary()
        # Define solution. If None, nothing will be added. Either dict or
        # regular solution
        if solution_fluxes is not None:
            if isinstance(solution_fluxes, Solution):
                solution_fluxes = solution_fluxes.fluxes.to_dict()
            json_dict.flux_solution = solution_fluxes
        # Get graph and add to json_dict
        json_dict.graph = self.graph.copy()
        json_dict.reaction_strings = {
            reaction: self.members.get_by_id(reaction).reaction
            for reaction in json_dict.graph.keys()
        }

        return json_dict.visualize(
            filepath=filename,
            vertical=self.vertical,
            color=[self.color_negative, self.color_positive],
            min_max=self.color_min_max,
            quantile=self.color_quantile,
            max_steps=self.color_max_steps,
            n_steps=self.color_n_steps,
        )


def model_convert(model: Model):
    """
    Converts the all Group objects in given model to proper cobramod
    :class:`cobramod.pathway.Pathway`
    """
    for index, group in enumerate(iterable=model.groups):
        model.groups[index] = Pathway._transform(obj=group)
