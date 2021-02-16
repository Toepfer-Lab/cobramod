"""Unittest for module pathway

This module includes the TestCase TestGroup. This checks the behaviour of the
new child of :obj:`cobra.core.group.Group` "Pathway". This new class is able
to use Escher for its visualizations.
"""
from json import loads
from logging import DEBUG
from pathlib import Path
from unittest import TestCase, main
from time import sleep

from cobra.core import DictList, Group

from cobramod.core import pathway as pt
from cobramod.core.extension import add_pathway
from cobramod.debug import debug_log
from cobramod.test import textbook_kegg, textbook, textbook_biocyc

# Debug must be set in level DEBUG for the test
debug_log.setLevel(DEBUG)
# Setting directory for data
dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class TestGroup(TestCase):
    def test___init__(self):
        # CASE 1: regular __init__
        test_group = pt.Pathway(id="test_group")
        self.assertIn(member="solution", container=dir(test_group))
        self.assertIn(member="members", container=dir(test_group))
        self.assertIn(member="kind", container=dir(test_group))
        self.assertIn(member="id", container=dir(test_group))
        self.assertIn(member="KIND_TYPES", container=dir(test_group))
        # CASE 2: two entities cannot be the same
        test_group2 = pt.Pathway(id="test_group")
        self.assertIsNot(expr1=test_group, expr2=test_group2)

    def test__filter(self):
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

    def test_solution(self):
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

    def test_add_members(self):
        test_model = textbook_kegg.copy()
        # CASE 1: Regular reactions
        test_group = pt.Pathway(id="test_group")
        test_group.add_members(
            new_members=[
                test_model.reactions.R00127_c,
                test_model.reactions.R00315_c,
            ]
        )
        test_list = [reaction.id for reaction in test_group.members]
        self.assertIn(member="R00127_c", container=test_list)
        self.assertIn(member="R00315_c", container=test_list)
        # CASE 2: Mixing reactions and metabolites
        test_group = pt.Pathway(id="test_group")
        self.assertRaises(
            TypeError,
            test_group.add_members,
            new_members=[
                test_model.reactions.R00127_c,
                test_model.metabolites.C00084_c,
            ],
        )

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
            test_group.add_members([reaction])
            test_model.add_groups([test_group])
        for index, item in enumerate(test_model.groups):
            test_model.groups[index] = pt.Pathway._transform(item)
        for group in test_model.groups:
            self.assertIsInstance(obj=group, cls=pt.Pathway)

    def test_visualize(self):
        # CASE 1: Members with initialization.
        test_model = textbook.copy()
        members = DictList()
        members.union(
            iterable=[
                test_model.reactions.EX_glc__D_e,
                test_model.reactions.GLCpts,
                test_model.reactions.G6PDH2r,
                test_model.reactions.PGL,
                test_model.reactions.GND,
            ]
        )
        test_group = pt.Pathway(id="test_group", members=members)
        test_group.graph = {
            "EX_glc__D_e": "GLCpts",
            "GLCpts": "G6PDH2r",
            "G6PDH2r": "PGL",
            "PGL": "GND",
            "GND": None,
        }
        test_builder = test_group.visualize()
        self.assertEqual(
            first=len(loads(test_builder.map_json)[1]["reactions"]), second=5
        )
        # CASE 2: Members after initialization.
        test_model = textbook.copy()
        test_group = pt.Pathway(id="test_group")
        for reaction in ("EX_glc__D_e", "GLCpts", "G6PDH2r", "PGL", "GND"):
            test_group.add_members(
                new_members=[test_model.reactions.get_by_id(reaction)]
            )
        self.assertEqual(
            first=len(loads(test_builder.map_json)[1]["reactions"]), second=5
        )
        # CASE 3: From copy of model.
        test_model.add_groups(group_list=[test_group])
        test_model_copy = test_model.copy()
        test_group = test_model_copy.groups.get_by_id("test_group")
        self.assertEqual(
            first=len(loads(test_builder.map_json)[1]["reactions"]), second=5
        )
        # CASE 4: Regular Biocyc
        test_model = textbook_biocyc.copy()
        add_pathway(
            model=test_model,
            pathway="SALVADEHYPOX-PWY",
            compartment="c",
            directory=dir_data,
            database="META",
            ignore_list=[],
            show_imbalance=False,
        )
        # Test fluxes
        test_pathway = test_model.groups.get_by_id("SALVADEHYPOX-PWY")
        self.assertEqual(first=len(test_pathway.members), second=5)
        test_solution = test_pathway.solution(solution=test_model.optimize())
        test_pathway.visualize(solution_fluxes=test_solution)
        sleep(1)
        add_pathway(
            model=test_model,
            pathway="PWY-1187",
            compartment="c",
            directory=dir_data,
            database="META",
            ignore_list=[],
            show_imbalance=False,
        )
        # Test fluxes
        test_pathway = test_model.groups.get_by_id("PWY-1187")
        self.assertEqual(first=len(test_pathway.members), second=14)
        test_solution = test_pathway.solution(solution=test_model.optimize())
        test_pathway.visualize(solution_fluxes=test_solution)


if __name__ == "__main__":
    main(verbosity=2)
