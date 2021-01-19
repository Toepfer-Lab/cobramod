"""Group-child class

The new class Pathway is child derived from :func:`cobra.core.group.Group`. It
extends some functionalities such as:

- solution: Obtain the solution for the specific members.
- visualize: get a :func:`escher.Builder` for that specific Pathway.
"""
from pathlib import Path
from typing import Dict, List

from escher import Builder
from cobra.core.group import Group
from cobra.core.solution import Solution
from cobra.core.reaction import Reaction
from cobra.core.dictlist import DictList
from pandas import Series

from cobramod.debug import debug_log
from cobramod.visualization.converter import JsonDictionary


class Pathway(Group):
    """
    A Sub-class from the original COBRApy :func:`cobra.Group`, which inherits
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
        parent :func:`cobra.core.group.Group`. Included is the order of the
        reactions.

        Args:
            id (str): Identifier of the pathway.
            name (str): Name for the patway.
            members (DictList): Includes the members for the pathway. It should
                include only :func`cobra.core.reaction.Reaction`. Defaults to
                None.
            kind (str): The kind of group, as specified for the Groups feature
                in the SBML level 3 package specification. Can be any of
                "classification", "partonomy", or "collection". The default is
                "collection". Please consult the SBML level 3 package
                specification to ensure you are using the proper value for
                kind. (Text extracted from parent.)
        """
        super().__init__(id=id, name=name, kind=kind, members=members)
        self.order: List[str] = list()
        # This loop has to be after __init__, otherwise, behaviour of members
        # changes
        for member in self.members:
            # Remove no Reactions
            if not isinstance(member, Reaction):
                self.remove_members(to_remove=[member])
                debug_log.debug(
                    f'Member "{member.id}" is not a Reaction. Skipped'
                )
            else:
                self.order.append(member.id)

    def _filter(self, solution: Solution, attribute: str) -> Series:
        """
        Filters given attribute of Solution. Only members of the Pathway will
        appear.
        """
        return getattr(solution, "fluxes").filter(
            items=[member.id for member in self.members], axis="index"
        )

    def add_members(self, new_members: List[Reaction]):
        """
        Add given list of :func:`cobra.core.reaction.Reaction` into the
        Pathway.

        Args:
            new_members:

        Raises:
            TypeError:
        """
        if not all([isinstance(member, Reaction) for member in new_members]):
            raise TypeError("Not all given members are Reactions.")
        super().add_members(new_members=new_members)
        # Extend order in order to use it later for the visualization.
        self.order.extend((reaction.id for reaction in new_members))

    def solution(self, solution: Solution) -> Solution:
        """

        Args:
            solution (Solution): Original COBRApy :func:`cobra.Solution` to
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
        Transform given :func:`cobra.Group` into a proper cobramod Pathway
        object.
        """
        return cls(
            id=obj.id, name=obj.name, members=obj.members, kind=obj.kind
        )

    def visualize(
        self,
        solution_fluxes: Dict[str, float] = None,
        filename: Path = None,
        canvas_height: float = 1500,
        canvas_width: float = 1500,
        **kwargs,
    ) -> Builder:
        """
        Returns a :func:`escher.Builder`, which can be use to create visual
        representations of the pathway.

        Args:
            solution_fluxes (dict): Dictionary with fluxes. The values will be
                then showed in the Builder. Defaults to None.
            filename (Path): Path for the HTML. If None is passed, then
                default to "pathway.html" in the current working directory.
            canvas_height (float): Height for the canvas. Defaults to 1500
            canvas_width (float): Width for the canvas. Defaults to 1500

        Keyword Argumengs:
            Read
            :func:`cobramod.visualization.converter.JsonDictionary.__init__`
        """
        # Set the canvas
        canvas = {
            "x": 0,
            "y": 0,
            "width": canvas_width,
            "height": canvas_height,
        }
        json_dict = JsonDictionary(canvas=canvas, **kwargs)
        member: str
        # This is to secure that the order of the members is respected, since
        # the order is important
        for member in self.order:
            reaction_string = self.members.get_by_id(member).reaction
            json_dict.add_reaction(string=reaction_string, identifier=member)
        # Define solution. If None, nothing will be added.
        if solution_fluxes is not None:
            json_dict.reaction_data = solution_fluxes
        return json_dict.visualize(filepath=filename)
