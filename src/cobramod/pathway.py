from cobra.core.group import Group
from cobra.core.solution import Solution
from cobra.core.dictlist import DictList
from pandas import Series


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
        members: DictList = None,
        kind: str = None,
    ):
        super().__init__(id=id, name=name, members=members, kind=kind)

    def _filter(self, solution: Solution, attribute: str) -> Series:
        """
        Filters given attribute of Solution. Only members of the Pathway will
        appear.
        """
        return getattr(solution, "fluxes").filter(
            items=[member.id for member in self.members], axis="index"
        )

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
        return cls(id=obj.id, name=obj.name, members=obj.name, kind=obj.kind)
