from contextlib import suppress
from typing import NamedTuple

from cobra.core.group import Group
from cobra.core.solution import Solution


class ReactionFlux(NamedTuple):
    reaction: str
    flux: float


class Pathway(Group):
    def __init__(self, id, name="", members=None, kind=None):
        super().__init__(id, name, members, kind)

    def fluxes(self):
        pass

    def _filter_fluxes(self, solution: Solution):
        for member in self.members:
            with suppress(KeyError):
                yield member.id, solution.fluxes[member.id]
