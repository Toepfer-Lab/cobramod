"""Unittest for module pathway

This module includes the TestCase TestGroup. This checks the behavior of the
new child of :obj:`cobra.core.group.Group` "Pathway". This new class is able
to use Escher for its visualizations.
"""
from json import loads
from logging import DEBUG
from pathlib import Path
from unittest import TestCase, main
from random import randint
from time import sleep

from cobra.core import DictList, Group
from cobra.io import read_sbml_model, write_sbml_model

from cobramod.core import pathway as pt
from cobramod.core.extension import add_pathway
from cobramod.core.pathway import Pathway
from cobramod.debug import debug_log
from cobramod.error import GraphKeyError
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
        self.assertDictEqual(
            d1=test_group.notes["ORDER"],
            d2={"R00315_c": None, "R00127_c": None},
        )
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
        for reaction in test_model.reactions[:5]:
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
        sleep(1)
        test_group.vertical = True
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
        # CASE 4a: Regular Biocyc
        test_model = textbook_biocyc.copy()
        add_pathway(
            model=test_model,
            pathway="SALVADEHYPOX-PWY",
            compartment="c",
            directory=dir_data,
            database="VCHO",
            ignore_list=[],
            show_imbalance=False,
        )
        # Test fluxes
        test_pathway = test_model.groups.get_by_id("SALVADEHYPOX-PWY")
        self.assertEqual(first=len(test_pathway.members), second=5)
        test_solution = {
            reaction.id: randint(-4, 4) for reaction in test_pathway.members
        }
        test_pathway.color_negative = "red"
        test_pathway.color_positive = "green"
        test_pathway.visualize(solution_fluxes=test_solution)
        sleep(1)
        test_pathway.vertical = True
        test_pathway.visualize(solution_fluxes=test_solution)
        sleep(1)
        # CASE 4b: Regular Biocyc
        add_pathway(
            model=test_model,
            pathway="PWY-1187",
            compartment="c",
            directory=dir_data,
            database="ARA",
            ignore_list=[],
            show_imbalance=False,
        )
        # Test fluxes
        test_pathway = test_model.groups.get_by_id("PWY-1187")
        self.assertEqual(first=len(test_pathway.members), second=14)
        self.assertTrue(expr=test_pathway.notes["ORDER"])
        test_solution = {
            reaction.id: randint(-4, 4) for reaction in test_pathway.members
        }
        test_pathway.color_negative = "purple"
        test_pathway.color_positive = "blue"
        test_pathway.visualize(solution_fluxes=test_solution)
        sleep(1)
        test_pathway.vertical = True
        test_pathway.visualize(solution_fluxes=test_solution)

    def test_model_convert(self):
        # CASE 1: regular conversion of Groups
        test_model = textbook.copy()
        for reaction in test_model.reactions[:5]:
            test_group = Group(id=reaction.id)
            test_model.add_groups(group_list=[test_group])
            test_model.groups.get_by_id(reaction.id).add_members(
                new_members=[reaction]
            )
        pt.model_convert(model=test_model)
        for group in test_model.groups:
            self.assertIsInstance(obj=group, cls=Pathway)
        # Test visualize
        test_model.groups[-1].visualize()
        sleep(1)

        # CASE 2: Regular Model
        filename = dir_input.joinpath("test_model02.sbml")
        test_model = read_sbml_model(str(filename))
        pt.model_convert(test_model)
        for group in test_model.groups:
            self.assertIsInstance(obj=group, cls=Pathway)
        test_model.groups[-1].visualize()
        sleep(1)

        # CASE 3: using add_pathway
        test_model = textbook_biocyc.copy()
        sequence = ["PEPDEPHOS-RXN", "PYRUVFORMLY-RXN", "FHLMULTI-RXN"]
        test_graph = {
            "PEPDEPHOS_RXN_c": "PYRUVFORMLY_RXN_c",
            "PYRUVFORMLY_RXN_c": "FHLMULTI_RXN_c",
            "FHLMULTI_RXN_c": None,
        }
        add_pathway(
            model=test_model,
            pathway=sequence,
            directory=dir_data,
            database="ECOLI",
            compartment="c",
            group="test_group",
        )

        filename = dir_input.joinpath("test_model.sbml")
        write_sbml_model(cobra_model=test_model, filename=str(filename))

        test_model = read_sbml_model(str(filename))
        pt.model_convert(test_model)
        test_pathway = test_model.groups.get_by_id("test_group")
        self.assertDictEqual(d1=test_graph, d2=test_pathway.notes["ORDER"])
        self.assertDictEqual(d1=test_graph, d2=test_pathway.graph)
        test_pathway.visualize()
        sleep(1)
        # CASE 4: Using a Group
        test_model = textbook_biocyc.copy()

        test_group = Group(id="curated_pathway")
        for reaction in ("GLCpts", "G6PDH2r", "PGL", "GND"):
            test_group.add_members([test_model.reactions.get_by_id(reaction)])
        test_model.add_groups([test_group])

        # Conversion to a Pathway
        pt.model_convert(model=test_model)
        test_model.groups.get_by_id("curated_pathway")
        test_model.groups.get_by_id("curated_pathway").vertical = True
        test_model.groups.get_by_id("curated_pathway").visualize()
        sleep(1)

    def test_modify_graph(self):
        test_model = textbook_kegg.copy()
        # CASE 1: merging test
        test_sequence = ["RXN-2206", "RXN-11414", "RXN-11422", "RXN-11430"]
        add_pathway(
            model=test_model,
            pathway=test_sequence,
            database="ARA",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
            avoid_list=["RXN-2206"],
        )
        test_sequence = ["RXN-11438", "RXN-2208", "RXN-2209", "RXN-2221"]
        add_pathway(
            model=test_model,
            pathway=test_sequence,
            database="ARA",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
        )
        test_group: Pathway = test_model.groups.get_by_id("custom_group")
        test_group.modify_graph(
            reaction="RXN_11430_c", next_reaction="RXN_11438_c"
        )
        self.assertDictContainsSubset(
            {"RXN_11430_c": "RXN_11438_c"}, test_group.graph
        )
        self.assertDictContainsSubset(
            {"RXN_11430_c": "RXN_11438_c"}, test_group.notes["ORDER"]
        )
        test_group.visualize()
        sleep(1)
        # CASE 2: Keys do not exists
        self.assertRaises(
            GraphKeyError,
            test_group.modify_graph,
            reaction="RXN_11438_c",
            next_reaction="NOT_RXN",
        )
        self.assertRaises(
            GraphKeyError,
            test_group.modify_graph,
            reaction="NOT_RXN",
            next_reaction="RXN_11438_c",
        )
        # CASE 3: Using None
        test_group.modify_graph(reaction="RXN_11430_c", next_reaction=None)
        self.assertDictContainsSubset({"RXN_11430_c": None}, test_group.graph)
        self.assertDictContainsSubset(
            {"RXN_11430_c": None}, test_group.notes["ORDER"]
        )


if __name__ == "__main__":
    main(verbosity=2)
