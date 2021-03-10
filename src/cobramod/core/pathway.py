#!/usr/bin/env python3
"""COBRApy Group-child class extension

The new class :class:`cobramod.pathway.Pathway" is child derived from
:class:`cobra.core.group.Group`. It extends some functionalities such as:

- solution: Obtain the solution for the specific members.
- visualize: get a :class:`escher.Builder` for that specific Pathway.
"""
from pathlib import Path
from typing import Dict, List, Union

from escher import Builder
from cobra.core.group import Group
from cobra.core.model import Model
from cobra.core.solution import Solution
from cobra.core.reaction import Reaction
from cobra.core.dictlist import DictList
from pandas import Series

from cobramod.debug import debug_log
from cobramod.visualization.converter import JsonDictionary


class Pathway(Group):
    """
    A Sub-class from the original COBRApy :class:`cobra.Group`, which inherits
    all its attributes and adds the method solution, to get a Solution for the
    members of this Class.
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
        # Loop has to be after __init__, otherwise, behaviour of class changes.
        # TODO: Is order necessary?
        self.order: List[str] = list()
        self.graph: dict = dict()
        for member in members:
            # Remove no Reactions
            if not isinstance(member, Reaction):
                debug_log.debug(
                    f'Member "{member.id}" is not a Reaction. Skipped'
                )
                continue
            self.add_members(new_members=[member])

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
        "order" is larger thatn zero. This check is made is to check if method
        :func:`cobra.core.model.Model.copy` is called. For copies, order must
        reset.
        """
        if len(self.members) == 0 and len(self.order) > 0:
            self.order: List[str] = list()
            self.graph: dict = dict()
            debug_log.debug(
                f'Attribute order from  pathway "{self.id}" reset. '
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

    def visualize(
        self,
        solution_fluxes: Union[Solution, Dict[str, float]] = None,
        filename: Path = None,
        vertical: bool = False,
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
        return json_dict.visualize(filepath=filename, vertical=vertical)


def model_convert(model: Model):
    """
    Converts the all Group objects in given model to proper cobramod
    :class:`cobramod.pathway.Pathway`
    """
    for index, group in enumerate(iterable=model.groups):
        model.groups[index] = Pathway._transform(obj=group)
