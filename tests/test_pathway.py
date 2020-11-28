import unittest

from cobra.core import DictList, Group

from cobramod import pathway as pt
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
                test_model.reactions.R00127_c,
                test_model.metabolites.C00084_c,
            ]
        )
        test_group = pt.Pathway(id="test_group", members=members)
        test_model.add_groups([test_group])
        # CASE 1: Fluxes
        test_serie = test_model.groups.get_by_id("test_group")._filter(
            solution=test_model.optimize(), attribute="fluxes"
        )
        self.assertEqual(first=test_serie["R00127_c"], second=0)
        # CASE 2: Shadows prices
        test_serie = test_model.groups.get_by_id("test_group")._filter(
            solution=test_model.optimize(), attribute="shadow_prices"
        )
        self.assertEqual(first=test_serie["R00127_c"], second=0)

    def test__solution(self):
        # Configuration
        test_model = textbook_kegg.copy()
        members = DictList()
        members.union(
            iterable=[
                test_model.reactions.R00127_c,
                test_model.reactions.R00315_c,
            ]
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
            self.assertTrue(expr="R00127_c" in attribute)
            self.assertTrue(expr="R00315_c" in attribute)

    def test__transform(self):
        # CASE 1: Outside the model
        test_model = textbook_kegg.copy()
        members = DictList()
        members.union(
            iterable=[
                test_model.reactions.R00127_c,
                test_model.reactions.R00315_c,
            ]
        )
        test_group = Group(id="test_group", members=members)
        test_group = pt.Pathway._transform(obj=test_group)
        self.assertIsInstance(obj=test_group, cls=pt.Pathway)
        # CASE 1a: Inside the model
        test_model = textbook_kegg.copy()
        members = DictList()
        members.union(
            iterable=[
                test_model.reactions.R00127_c,
                test_model.reactions.R00315_c,
            ]
        )
        test_group = Group(id="test_group", members=members)
        test_model.add_groups([test_group])
        test_model.groups[0] = pt.Pathway._transform(
            obj=test_model.groups.get_by_id("test_group")
        )
        self.assertIsInstance(
            obj=test_model.groups.get_by_id("test_group"), cls=pt.Pathway
        )
        # CASE 1b: Inside the model, multiple
        test_model = textbook_kegg.copy()
        for reaction in test_model.reactions:
            test_group = Group(id=reaction.id)
            test_model.add_groups([test_group])
        for index, item in enumerate(test_model.groups):
            test_model.groups[index] = pt.Pathway._transform(item)
        for group in test_model.groups:
            self.assertIsInstance(obj=group, cls=pt.Pathway)


if __name__ == "__main__":
    unittest.main(verbosity=2)
