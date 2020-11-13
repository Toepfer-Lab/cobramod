import unittest

from cobra.core.dictlist import DictList

import cobramod.pathway as pt
from cobramod.test import textbook_kegg


class TestGroup(unittest.TestCase):
    def test___init__(self):
        # CASE 1: regular __init__
        test_group = pt.Pathway(id="test_group")
        self.assertIn(member="solution", container=dir(test_group))
        self.assertIn(member="members", container=dir(test_group))
        self.assertIn(member="kind", container=dir(test_group))
        self.assertIn(member="id", container=dir(test_group))
        self.assertIn(member="KIND_TYPES", container=dir(test_group))

    def test__filter(self):
        # Configuration
        test_model = textbook_kegg.copy()
        members = DictList()
        members.union(
            iterable=[
                test_model.reactions.ADK1,
                test_model.metabolites.C00084_c,
            ]
        )
        test_group = pt.Pathway(id="test_group", members=members)
        test_model.add_groups([test_group])
        # CASE 1: Fluxes
        test_serie = test_model.groups.get_by_id("test_group")._filter(
            solution=test_model.optimize(), attribute="fluxes"
        )
        self.assertEqual(first=test_serie["ADK1"], second=0)
        # CASE 2: Shadows prices
        test_serie = test_model.groups.get_by_id("test_group")._filter(
            solution=test_model.optimize(), attribute="shadow_prices"
        )
        self.assertEqual(first=test_serie["ADK1"], second=0)

    def test__solution(self):
        # Configuration
        test_model = textbook_kegg.copy()
        members = DictList()
        members.union(
            iterable=[test_model.reactions.ADK1, test_model.reactions.ACKr]
        )
        test_group = pt.Pathway(id="test_group", members=members)
        test_model.add_groups([test_group])

        test_solution = test_model.groups.get_by_id("test_group").solution(
            solution=test_model.optimize()
        )
        # CASE 1: Regular fluxes
        for attribute in (
            test_solution.fluxes,
            test_solution.reduced_costs,
            test_solution.shadow_prices,
        ):
            self.assertTrue(expr="ADK1" in attribute)
            self.assertTrue(expr="ACKr" in attribute)


if __name__ == "__main__":
    unittest.main(verbosity=2)
