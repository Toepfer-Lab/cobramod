import unittest

from cobra.core.dictlist import DictList

import cobramod.pathway as pt
from cobramod.test import textbook_kegg


class TestGroup(unittest.TestCase):
    def test__filter_fluxes(self):
        # CASE 1: Regular Group
        test_model = textbook_kegg.copy()
        # test_model.groups = gr.Pathway
        members = DictList()
        members.union(
            iterable=[
                test_model.reactions.ADK1,
                test_model.metabolites.C00084_c,
            ]
        )
        test_group = pt.Pathway(id="test_group", members=members)
        test_model.add_groups([test_group])
        test_solution = test_model.optimize()
        list(
            test_model.groups.get_by_id("test_group")._filter_fluxes(
                solution=test_solution
            )
        )

        # gr._filter_fluxes(
        #     solution=test_solution,
        #     group=test_model.groups.get_by_id("test_group"),
        # )
        pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
