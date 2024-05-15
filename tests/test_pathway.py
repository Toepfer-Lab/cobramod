#!/usr/bin/env python3
"""Unittest for module pathway

This module includes the TestCase TestGroup. This checks the behavior of the
new child of :obj:`cobra.core.group.Group` "Pathway". This new class is able
to use Escher for its visualizations.
"""

import unittest
from json import loads
from pathlib import Path
from random import randint
from tempfile import TemporaryDirectory
from time import sleep
from webbrowser import open as web_open

import cobra.core as cobra_core
import cobra.io as cobra_io
import cobramod.error as cmod_error
import cobramod.test as cmod_test
import pandas as pd
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.core import pathway as pt
from cobramod.core.extension import add_pathway
from cobramod.debug import change_to_debug

change_to_debug()

dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class TestGroup(unittest.TestCase):
    def test___init__(self):
        # CASE: regular __init__
        test_group = pt.Pathway(id="test_group")
        self.assertIn(member="solution", container=dir(test_group))
        self.assertIn(member="members", container=dir(test_group))
        self.assertIn(member="kind", container=dir(test_group))
        self.assertIn(member="id", container=dir(test_group))
        self.assertIn(member="KIND_TYPES", container=dir(test_group))

        # CASE: two entities cannot be the same
        test_group2 = pt.Pathway(id="test_group")
        self.assertIsNot(expr1=test_group, expr2=test_group2)

    def test__filter(self):
        # Configuration
        test_model = cmod_test.textbook_kegg.copy()
        members = {
            test_model.reactions.R00127_c,
            test_model.reactions.R00315_c,
        }
        test_group = pt.Pathway(id="test_group", members=members)
        test_model.add_groups([test_group])

        # CASE: Fluxes
        test_pathway = test_model.groups.get_by_id("test_group")
        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        test_serie = test_pathway._filter(
            solution=test_model.optimize(), attribute="fluxes"
        )
        self.assertEqual(first=test_serie["R00127_c"], second=0)

        # CASE: reduced prices
        test_pathway = test_model.groups.get_by_id("test_group")
        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        test_serie = test_pathway._filter(
            solution=test_model.optimize(), attribute="reduced_costs"
        )
        self.assertEqual(first=test_serie["R00127_c"], second=0)

    def test_solution(self):
        # Configuration
        test_model = cmod_test.textbook_kegg.copy()
        members = {
            test_model.reactions.R00127_c,
            test_model.reactions.R00315_c,
        }

        test_group = pt.Pathway(id="test_group", members=members)
        test_model.add_groups([test_group])
        test_solution = test_group.solution(solution=test_model.optimize())
        # CASE: Regular fluxes

        for attribute in (
            test_solution.fluxes,
            test_solution.reduced_costs,
        ):
            if not isinstance(attribute, pd.Series):
                raise TypeError("Not a valid pandas.Series")

            self.assertTrue(expr="R00127_c" in attribute)
            self.assertTrue(expr="R00315_c" in attribute)

    def test_add_members(self):
        test_model = cmod_test.textbook_kegg.copy()
        # CASE: Regular reactions
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

        # CASE: Mixing reactions and metabolites
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
        # CASE: Outside the model
        test_model = cmod_test.textbook_kegg.copy()
        members = {
            test_model.reactions.R00127_c,
            test_model.reactions.R00315_c,
        }
        test_group = cobra_core.Group(id="test_group", members=members)
        test_group = pt.Pathway._transform(obj=test_group)
        self.assertIsInstance(obj=test_group, cls=pt.Pathway)

        # CASE: Inside the model
        test_model = cmod_test.textbook_kegg.copy()
        members = {
            test_model.reactions.R00127_c,
            test_model.reactions.R00315_c,
        }
        test_group = cobra_core.Group(id="test_group", members=members)
        test_model.add_groups([test_group])
        test_model.groups[0] = pt.Pathway._transform(obj=test_group)

        self.assertIsInstance(obj=test_model.groups[0], cls=pt.Pathway)

        # CASE: Inside the model, multiple
        test_model = cmod_test.textbook_kegg.copy()

        test_iterable = [
            test_model.reactions[0],
            test_model.reactions[1],
            test_model.reactions[2],
            test_model.reactions[3],
            test_model.reactions[4],
        ]
        for reaction in test_iterable:
            test_group = cobra_core.Group(id=reaction.id)
            test_group.add_members([reaction])
            test_model.add_groups([test_group])

        for index, item in enumerate(test_model.groups):
            test_model.groups[index] = pt.Pathway._transform(item)

        for group in test_model.groups:
            self.assertIsInstance(obj=group, cls=pt.Pathway)

    def test_visualize(self):
        # CASE: Members with initialization.
        test_model = cmod_test.textbook.copy()
        members = {
            test_model.reactions.EX_glc__D_e,
            test_model.reactions.GLCpts,
            test_model.reactions.G6PDH2r,
            test_model.reactions.PGL,
            test_model.reactions.GND,
        }
        test_group = pt.Pathway(id="test_group", members=members)
        test_group.graph = {
            "EX_glc__D_e": "GLCpts",
            "GLCpts": "G6PDH2r",
            "G6PDH2r": "PGL",
            "PGL": "GND",
            "GND": None,
        }
        test_builder = test_group.visualize(vis="escher-custom")

        web_open("pathway.html")
        sleep(1)
        test_group.vertical = True
        test_builder = test_group.visualize(vis="escher-custom")

        self.assertEqual(
            first=len(loads(test_builder.map_json)[1]["reactions"]),  # type: ignore
            second=5,
        )
        # CASE: Members after initialization.
        test_model = cmod_test.textbook.copy()
        test_group = pt.Pathway(id="test_group")
        for identifier in ("EX_glc__D_e", "GLCpts", "G6PDH2r", "PGL", "GND"):
            reaction = test_model.reactions.get_by_id(identifier)
            if not isinstance(reaction, cobra_core.Reaction):
                raise TypeError("Not a valid COBRApy object")

            test_group.add_members(new_members=[reaction])
        self.assertEqual(
            first=len(loads(test_builder.map_json)[1]["reactions"]),  # type: ignore
            second=5,
        )
        # CASE: From copy of model.
        test_model.add_groups(group_list=[test_group])
        test_model_copy = test_model.copy()
        test_group = test_model_copy.groups.get_by_id("test_group")
        self.assertEqual(
            first=len(loads(test_builder.map_json)[1]["reactions"]),  # type: ignore
            second=5,
        )
        # CASE: Regular Biocyc
        test_model = cmod_test.textbook_biocyc.copy()
        add_pathway(
            model=test_model,
            pathway="SALVADEHYPOX-PWY",
            compartment="c",
            directory=dir_data,
            database="ECO",
            ignore_list=[],
            show_imbalance=False,
        )
        # Test fluxes
        test_pathway = test_model.groups.get_by_id("SALVADEHYPOX-PWY")
        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        self.assertEqual(first=len(test_pathway.members), second=5)
        test_solution = {
            reaction.id: randint(-4, 4) for reaction in test_pathway.members
        }
        test_pathway.color_negative = "red"
        test_pathway.color_positive = "green"
        test_pathway.visualize(
            solution_fluxes=test_solution, vis="escher-custom"
        )
        web_open("pathway.html")
        sleep(1)

        test_pathway.vertical = True
        test_pathway.visualize(
            solution_fluxes=test_solution, vis="escher-custom"
        )
        web_open("pathway.html")
        sleep(1)

        # CASE: Regular Biocyc
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
        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        self.assertEqual(first=len(test_pathway.members), second=14)
        self.assertTrue(expr=test_pathway.notes["ORDER"])
        test_solution: dict[str, float] = {
            reaction.id: randint(-4, 4) for reaction in test_pathway.members
        }
        test_pathway.color_negative = "purple"
        test_pathway.color_positive = "blue"
        test_pathway.visualize(
            solution_fluxes=test_solution, vis="escher-custom"
        )
        web_open("pathway.html")
        sleep(1)
        test_pathway.vertical = True
        test_pathway.visualize(
            solution_fluxes=test_solution, vis="escher-custom"
        )

    def test_model_convert(self):
        # CASE: regular conversion of Groups
        test_model = cmod_test.textbook.copy()
        test_list = [
            test_model.reactions[0],
            test_model.reactions[1],
            test_model.reactions[2],
            test_model.reactions[3],
            test_model.reactions[4],
        ]
        for reaction in test_list:
            test_group = cobra_core.Group(id=reaction.id)
            test_model.add_groups(group_list=[test_group])

            test_group.add_members(new_members=[reaction])

        pt.model_convert(model=test_model)
        for group in test_model.groups:
            self.assertIsInstance(obj=group, cls=pt.Pathway)

        # Test visualize
        test_pathway = test_model.groups[-1]
        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")
        test_pathway.visualize()
        web_open("pathway.html")
        sleep(1)

        # CASE: Regular Model
        filename = dir_input.joinpath("test_model.sbml")
        test_model = cobra_io.read_sbml_model(str(filename))
        pt.model_convert(test_model)
        for group in test_model.groups:
            self.assertIsInstance(obj=group, cls=pt.Pathway)
        test_pathway = test_model.groups[-1]

        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        test_pathway.visualize()
        web_open("pathway.html")
        sleep(1)

        # CASE: using add_pathway
        test_model = cmod_test.textbook_biocyc.copy()
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

        with TemporaryDirectory() as tmp_dir:
            filename = Path(tmp_dir).joinpath("test_model.sbml")
            cobra_io.write_sbml_model(
                cobra_model=test_model, filename=str(filename)
            )

            test_model = cobra_io.read_sbml_model(str(filename))

        pt.model_convert(test_model)
        test_pathway = test_model.groups.get_by_id("test_group")

        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        self.assertDictEqual(d1=test_graph, d2=test_pathway.notes["ORDER"])
        self.assertDictEqual(d1=test_graph, d2=test_pathway.graph)
        test_pathway.visualize()
        web_open("pathway.html")
        sleep(1)

        # CASE: Using a Group
        test_model = cmod_test.textbook_biocyc.copy()

        test_group = cobra_core.Group(id="curated_pathway")
        for reaction in ("GLCpts", "G6PDH2r", "PGL", "GND"):
            test_group.add_members([test_model.reactions.get_by_id(reaction)])
        test_model.add_groups([test_group])

        # Conversion to a Pathway
        pt.model_convert(model=test_model)
        test_pathway = test_model.groups.get_by_id("curated_pathway")
        if not isinstance(test_pathway, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        test_pathway.vertical = True
        test_pathway.visualize()
        web_open("pathway.html")

    def test_modify_graph(self):
        test_model = cmod_test.textbook_kegg.copy()
        # CASE: merging test
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
        test_group = test_model.groups.get_by_id("custom_group")

        if not isinstance(test_group, pt.Pathway):
            raise TypeError("Not a valid CobraMod Pathway")

        test_group.modify_graph(
            reaction="RXN_11430_c", next_reaction="RXN_11438_c"
        )
        subset = {"RXN_11430_c": "RXN_11438_c"}
        self.assertEqual(test_group.graph, test_group.graph | subset)
        self.assertEqual(test_group.graph, test_group.notes["ORDER"] | subset)
        test_group.visualize()
        web_open("pathway.html")
        sleep(1)

        # CASE: Keys do not exists
        self.assertRaises(
            cmod_error.GraphKeyError,
            test_group.modify_graph,
            reaction="RXN_11438_c",
            next_reaction="NOT_RXN",
        )
        self.assertRaises(
            cmod_error.GraphKeyError,
            test_group.modify_graph,
            reaction="NOT_RXN",
            next_reaction="RXN_11438_c",
        )
        # CASE: Using None

        subset = {"RXN_11430_c": None}
        test_group.modify_graph(reaction="RXN_11430_c", next_reaction=None)
        self.assertEqual(test_group.graph, test_group.graph | subset)
        self.assertEqual(
            test_group.notes["ORDER"], test_group.notes["ORDER"] | subset
        )


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
