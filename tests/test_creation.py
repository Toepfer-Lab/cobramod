#!/usr/bin/env python3
"""Unittest for module creation

In this test, multiple single functions for the creation of reactions and
metabolites are tested. Important is the independency of the functions, unless
they are part of the API such as add_metabolites. Therefore the modules is
separated in:

- SimpleFunctions: Creation of objects
- ComplexFunctions: Functions, that uses multiple simple functions.
"""

import unittest
from logging import DEBUG
from pathlib import Path

import cobra.core as cobra_core
import cobramod.retrieval as cmod_retrieval
import requests
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.core import creation as cr
from cobramod.debug import debug_log
from cobramod.error import NoDelimiter, WrongSyntax
from cobramod.parsing.db_version import DataVersionConfigurator
from cobramod.test import textbook_kegg

NAME = "test_model"

debug_log.setLevel(DEBUG)

dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class SimpleFunctions(unittest.TestCase):
    """
    Test for simple test such as creating metabolites from string or from
    files.
    """

    @classmethod
    def setUp(cls):
        # ToDo Versioning should be tested in its own test class
        data_conf = DataVersionConfigurator()
        data_conf.ignore_db_versions = True

    def test_find_replacements(self):
        replace_dict = {
            "ACETALD-DEHYDROG-RXN": "R00228_c",
            "WATER_c": "C00001_c",
            "WATER_d": "C00001_d",
        }
        # CASE 1: Reaction
        test_replacement = cr.find_replacements(
            identifier="WATER_c",
            obj_type="metabolites",
            model=textbook_kegg,
            replace_dict=replace_dict,
        )
        self.assertIsInstance(obj=test_replacement, cls=cobra_core.Metabolite)
        # CASE: String
        test_replacement = cr.find_replacements(
            identifier="WATER_d",
            obj_type="metabolites",
            model=textbook_kegg,
            replace_dict=replace_dict,
        )
        self.assertIsNone(test_replacement)
        # CASE : String
        test_replacement = cr.find_replacements(
            identifier="ACETALD-DEHYDROG-RXN",
            obj_type="reactions",
            model=textbook_kegg,
            replace_dict=replace_dict,
        )
        self.assertIsInstance(obj=test_replacement, cls=cobra_core.Reaction)

    def test_metabolite_from_string(self):
        # CASE 1: Correct input
        test_string = "MALTOSE_b, MALTOSE[b], b, C12H22O11, 0"
        test_metabolite = cr.metabolite_from_string(
            line_string=test_string,
            model=cobra_core.Model(NAME),
            replacement={},
        )
        # Checking that new meta is in Model
        self.assertIsInstance(obj=test_metabolite, cls=cobra_core.Metabolite)
        self.assertEqual(first=test_metabolite.id, second="MALTOSE_b")
        self.assertEqual(first=test_metabolite.name, second="MALTOSE[b]")
        self.assertEqual(first=test_metabolite.charge, second=0)
        self.assertEqual(first=test_metabolite.formula, second="C12H22O11")

    def test_metabolite_from_data(self):
        # CASE: regular META
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="HOMOMETHIONINE",
            database="META",
            # debug_level=10,
        )
        test_metabolite = cr.metabolite_from_data(data, compartment="c")

        self.assertIsInstance(obj=test_metabolite, cls=cobra_core.Metabolite)
        self.assertEqual(first=test_metabolite.id, second="HOMOMETHIONINE_c")
        self.assertEqual(first=test_metabolite.name, second="L-homomethionine")
        self.assertEqual(first=test_metabolite.formula, second="C6H13N1O2S1")
        # CASE: Find translation
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="WATER",
            database="META",
            # debug_level=10,
        )
        test_metabolite = cr.metabolite_from_data(
            data, compartment="c", model=textbook_kegg
        )

        self.assertEqual(first=test_metabolite.id, second="C00001_c")
        # CASE: Find translation (BIGG)
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="h2o",
            database="BIGG",
            # debug_level=10,
            model_id="universal",
        )
        test_metabolite = cr.metabolite_from_data(
            data, compartment="c", model=textbook_kegg
        )
        self.assertEqual(first=test_metabolite.id, second="C00001_c")

    def test_convert_string_metabolite(self):
        # CASE 1: retrieval from META
        test_metabolite = cr.convert_string_metabolite(
            line="HOMOMETHIONINE, c",
            model=cobra_core.Model(NAME),
            directory=dir_data,
            database="META",
            replacement={},
            model_id=None,
        )
        self.assertIsInstance(obj=test_metabolite, cls=cobra_core.Metabolite)
        self.assertEqual(first=test_metabolite.id, second="HOMOMETHIONINE_c")
        self.assertEqual(first=test_metabolite.name, second="L-homomethionine")
        self.assertEqual(first=test_metabolite.formula, second="C6H13N1O2S1")
        self.assertEqual(first=test_metabolite.compartment, second="c")
        # CASE: custom metabolite
        test_metabolite = cr.convert_string_metabolite(
            line="MALTOSE_b, MALTOSE[b], b, C12H22O11, 0",
            model=cobra_core.Model(NAME),
            replacement={},
            directory=dir_data,
            database=None,
            model_id=None,
        )
        self.assertEqual(first=test_metabolite.id, second="MALTOSE_b")
        self.assertEqual(first=test_metabolite.name, second="MALTOSE[b]")
        self.assertEqual(first=test_metabolite.charge, second=0)
        self.assertEqual(first=test_metabolite.formula, second="C12H22O11")
        # CASE: custom metabolite with replacement
        test_metabolite = cr.convert_string_metabolite(
            line="MALTOSE_b, MALTOSE[c], b, C12H22O11, 0",
            model=cobra_core.Model(NAME),
            replacement={"MALTOSE_b": "MALTOSE_c"},
            directory=dir_data,
            database=None,
            model_id=None,
        )
        self.assertEqual(first=test_metabolite.id, second="MALTOSE_c")
        self.assertEqual(first=test_metabolite.name, second="MALTOSE[c]")
        self.assertEqual(first=test_metabolite.charge, second=0)
        self.assertEqual(first=test_metabolite.formula, second="C12H22O11")

    def test_get_file_metabolites(self):
        test_model = cobra_core.Model(NAME)
        # CASE: Testing if file is not found
        self.assertRaises(
            FileNotFoundError,
            cr.get_file_metabolites,
            test_model,
            Path.cwd().joinpath("no_file"),
            replacement={},
            directory=dir_data,
            database=None,
            model_id=None,
        )
        # CASE: Metabolite is not found (or misspelled)
        self.assertRaises(
            requests.HTTPError,
            cr.get_file_metabolites,
            test_model,
            dir_input.joinpath("metabolites_02_misspelled.txt"),
            replacement={},
            directory=dir_data,
            database="META",
            model_id=None,
        )
        # CASE: Bad format (e. g. charge is missing).
        self.assertRaises(
            WrongSyntax,
            cr.get_file_metabolites,
            model=test_model,
            filename=dir_input.joinpath("metabolites_03_badFormat.txt"),
            directory=dir_data,
            database="META",
            replacement={},
            model_id=None,
        )
        # CASE: Normal input
        test_list = cr.get_file_metabolites(
            model=test_model,
            filename=dir_input.joinpath("metabolites_01_normal.txt"),
            directory=dir_data,
            database="META",
            replacement={},
            model_id=None,
        )
        self.assertEqual(len(test_list), 2)

        # Both Names in metabolites
        test_names = ["HOMOMETHIONINE_c", "MALTOSE_b"]
        self.assertListEqual(
            list1=[member.id for member in test_list], list2=test_names
        )

    def test_get_reaction(self):
        # CASE: Regular Biocyc reaction
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="OXALODECARB-RXN",
            database="ARA",
            # debug_level=10,
        )
        test_reaction = cr.get_reaction(
            data=data,
            compartment="c",
            replacement={},
            model=cobra_core.Model(NAME),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIsInstance(obj=test_reaction, cls=cobra_core.Reaction)
        self.assertTupleEqual(tuple1=test_reaction.bounds, tuple2=(0, 1000))
        test_list = [
            "PYRUVATE_c",
            "CARBON_DIOXIDE_c",
            "OXALACETIC_ACID_c",
            "PROTON_c",
        ]
        self.assertCountEqual(
            [meta.id for meta in test_reaction.metabolites],
            test_list,
        )
        # CASE: Regular KEGG reaction in compartment "p". Irreversible
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="R02736",
            database="KEGG",
            # debug_level=10,
            # genome="hsa",
        )
        test_reaction = cr.get_reaction(
            data=data,
            compartment="p",
            replacement={},
            model=cobra_core.Model(NAME),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIsInstance(obj=test_reaction, cls=cobra_core.Reaction)
        # NOTE: dropping support from comments
        self.assertTupleEqual(tuple1=test_reaction.bounds, tuple2=(-1000, 1000))
        test_list = [
            "C01172_p",
            "C00006_p",
            "C01236_p",
            "C00005_p",
            "C00080_p",
        ]
        self.assertCountEqual(
            [meta.id for meta in test_reaction.metabolites],
            test_list,
        )
        # CASE: Reactions with different coefficients
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="R00114",
            database="KEGG",
            # debug_level=10,
            # genome="eco",
        )
        test_reaction = cr.get_reaction(
            data=data,
            compartment="p",
            replacement={},
            model=cobra_core.Model(NAME),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIsInstance(obj=test_reaction, cls=cobra_core.Reaction)
        self.assertTupleEqual(tuple1=test_reaction.bounds, tuple2=(-1000, 1000))
        test_list = [
            "C00025_p",
            "C00006_p",
            "C00064_p",
            "C00026_p",
            "C00080_p",
            "C00005_p",
        ]
        self.assertListEqual(
            list1=sorted([meta.id for meta in test_reaction.metabolites]),
            list2=sorted(test_list),
        )
        self.assertEqual(
            first=test_reaction.get_coefficient("C00025_p"), second=-2
        )
        # CASE: transport reaction, simple
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="TRANS-RXN0-574",
            database="META",
            # debug_level=10,
        )
        test_reaction = cr.get_reaction(
            data=data,
            compartment="p",
            replacement={},
            model=cobra_core.Model(NAME),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIn(
            member="Glucopyranose_p",
            container=[item.id for item in test_reaction.metabolites],
        )
        self.assertIn(
            member="Glucopyranose_e",
            container=[item.id for item in test_reaction.metabolites],
        )
        # CASE: transport reaction, complex
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="ABC-56-RXN",
            database="META",
            # debug_level=10,
        )
        test_reaction = cr.get_reaction(
            data=data,
            compartment="p",
            replacement={},
            model=cobra_core.Model(NAME),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIn(
            member="Aliphatic_Sulfonates_e",
            container=[item.id for item in test_reaction.metabolites],
        )
        self.assertIn(
            member="Aliphatic_Sulfonates_p",
            container=[item.id for item in test_reaction.metabolites],
        )
        self.assertEqual(
            first=test_reaction.get_coefficient("Aliphatic_Sulfonates_e"),
            second=-1,
        )
        # CASE: Retrieve reaction with translated equivalents
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="ACALD",
            database="BIGG",
            # debug_level=10,
            model_id="e_coli_core",
        )
        test_reaction = cr.get_reaction(
            data=data,
            compartment="c",
            replacement={},
            model=textbook_kegg,
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertEqual(first=test_reaction.id, second="R00228_c")
        # CASE: Retrieve reaction with translated equivalents. Check for
        # metabolites.
        data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="ADENODEAMIN-RXN",
            database="META",
            # debug_level=10,
        )
        test_reaction = cr.get_reaction(
            data=data,
            compartment="c",
            replacement={},
            model=textbook_kegg,
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertEqual(first=test_reaction.id, second="ADENODEAMIN_RXN_c")
        # WATER
        self.assertIn(
            member="C00001_c",
            container=[
                metabolite.id for metabolite in test_reaction.metabolites
            ],
        )
        # PROTON
        self.assertIn(
            member="C00080_c",
            container=[
                metabolite.id for metabolite in test_reaction.metabolites
            ],
        )

    def test_convert_string_reaction(self):
        # CASE: using delimiter, compartment is cytosol
        test_model = cobra_core.Model(NAME)
        test_line = (
            "RXN_17742_c, RXN_17742_c | "
            "1 Oxidized-ferredoxins_c <-> 1 Reduced-ferredoxins_c "
        )
        test_reaction = cr.string_to_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
            replacement={},
            model_id=None,
            genome=None,
        )
        self.assertEqual(first="RXN_17742_c", second=test_reaction.id)
        for metabolite in ["Oxidized_ferredoxins_c", "Reduced_ferredoxins_c"]:
            self.assertIn(
                member=metabolite,
                container=[meta.id for meta in test_reaction.metabolites],
            )
        # CASE: No delimiter
        test_model = cobra_core.Model(NAME)
        test_line = "RXN-14462, c"
        test_reaction = cr.string_to_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
            replacement={},
            model_id=None,
            genome=None,
        )
        self.assertEqual(first="RXN_14462_c", second=test_reaction.id)

        # CASE: No delimiter, compartment p
        test_model = cobra_core.Model(NAME)
        test_line = "RXN-14462, p"
        test_reaction = cr.string_to_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="GCF_000020025",
            stop_imbalance=False,
            show_imbalance=True,
            replacement={},
            model_id=None,
            genome=None,
        )
        self.assertEqual(first="RXN_14462_p", second=test_reaction.id)

        # CASE: 100% Custom metabolite
        test_model = cobra_core.Model(NAME)
        test_line = "CUSTOM_rxn1_p, Custom_reaction | Meta_A_p --> Meta_B_p"
        test_reaction = cr.string_to_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
            replacement={},
            model_id=None,
            genome=None,
        )
        self.assertEqual(first="CUSTOM_rxn1_p", second=test_reaction.id)
        for metabolite in ["Meta_A_p", "Meta_B_p"]:
            self.assertIn(
                member=metabolite,
                container=[meta.id for meta in test_reaction.metabolites],
            )

    def test_reaction_from_string(self):
        # CASE: wrong format, no delimiter
        self.assertRaises(
            NoDelimiter,
            cr.reaction_from_string,
            line_string="GLC_cb, GLC_cb GLC_c <-> GLC_b",
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
            model_id=None,
            replacement={},
        )

        # CASE: Normal, ID and name differ
        test_reaction = cr.reaction_from_string(
            line_string="GLC_cb, Glucose Transport| 2 GLC_c <=> GLC_b",
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=False,
            replacement={},
            model_id=None,
        )
        self.assertEqual(first=test_reaction.id, second="GLC_cb")
        self.assertEqual(first=test_reaction.name, second="Glucose Transport")
        self.assertIsInstance(obj=test_reaction, cls=cobra_core.Reaction)

        test_list = ["GLC_c", "GLC_b"]
        for metabolite in test_list:
            self.assertIn(
                member=metabolite,
                container=[meta.id for meta in test_reaction.metabolites],
            )
        self.assertEqual(
            first=-2, second=test_reaction.get_coefficient("GLC_c")
        )
        self.assertEqual(first=1, second=test_reaction.get_coefficient("GLC_b"))

        # CASE: Larger reaction
        test_reaction = cr.reaction_from_string(
            line_string="GLC_cb, Glucose Transport|"
            + "4 GLT_c + 2 GLC_c <=> GLC_b + 2 GLT_b",
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=False,
            replacement={},
            model_id=None,
        )
        self.assertEqual(
            first=-4, second=test_reaction.get_coefficient("GLT_c")
        )
        self.assertEqual(first=2, second=test_reaction.get_coefficient("GLT_b"))

    def test_get_file_reactions(self):
        test_model = cobra_core.Model(NAME)
        test_list = cr.get_file_reactions(
            model=test_model,
            filename=dir_input.joinpath("reactions_normal.txt"),
            directory=dir_data,
            database="META",
            replacement={},
            model_id=None,
            genome=None,
            stop_imbalance=False,
            show_imbalance=True,
        )
        test_names = ["GLC_cb", "RXN_14462_p"]
        self.assertListEqual(
            list1=test_names, list2=[reaction.id for reaction in test_list]
        )


class ComplexFunctions(unittest.TestCase):
    def test_create_object(self):
        # CASE: metabolite from MetaCyc
        test_object = cr.create_object(
            identifier="GLC",
            directory=dir_data,
            compartment="c",
            database="META",
        )
        self.assertIsInstance(obj=test_object, cls=cobra_core.Metabolite)
        self.assertEqual(first=getattr(test_object, "compartment"), second="c")
        # CASE: metabolite from KEGG
        test_object = cr.create_object(
            identifier="C00026",
            directory=dir_data,
            compartment="p",
            database="KEGG",
        )
        self.assertIsInstance(obj=test_object, cls=cobra_core.Metabolite)
        self.assertEqual(first=getattr(test_object, "compartment"), second="p")
        # CASE: reaction from MetaCyc
        test_object = cr.create_object(
            identifier="2.6.1.58-RXN",
            directory=dir_data,
            compartment="p",
            database="META",
        )
        self.assertIsInstance(obj=test_object, cls=cobra_core.Reaction)
        self.assertIn(
            member="p", container=getattr(test_object, "compartments")
        )
        # CASE: reaction from KEGG
        test_object = cr.create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            genome="hsa",
        )
        self.assertIsInstance(obj=test_object, cls=cobra_core.Reaction)
        self.assertIn(
            member="c", container=getattr(test_object, "compartments")
        )
        self.assertCountEqual(
            first=[gene.id for gene in getattr(test_object, "genes")],
            second=["9563", "2539"],
        )

        # CASE: Reaction of Biocyc which could be a pathway as well
        test_object = cr.create_object(
            identifier="AMONITRO-RXN",
            directory=dir_data,
            compartment="c",
            database="META",
            show_imbalance=False,
        )
        self.assertIsInstance(obj=test_object, cls=cobra_core.Reaction)
        # CASE: Metabolites has correct formula
        test_object = cr.create_object(
            identifier="C00404",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            show_imbalance=False,
        )
        if not isinstance(test_object, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy object")
        if test_object.formula is None:
            raise TypeError("Formula is empty!")

        self.assertEqual(test_object.formula.find("("), -1)

        # CASE: Replacing names (Metabolite)
        test_replacement = {
            "GLC": "GLUCOSE",
            "AMONITRO-RXN": "monooxygenase",
            "Donor-H2": "Donor_new_name",
        }

        test_object = cr.create_object(
            identifier="GLC",
            directory=dir_data,
            compartment="c",
            database="META",
            replacement=test_replacement,
        )
        self.assertIsInstance(obj=test_object, cls=cobra_core.Metabolite)
        self.assertEqual(first=getattr(test_object, "compartment"), second="c")
        self.assertEqual(first=getattr(test_object, "id"), second="GLUCOSE_c")

        # CASE: Replacing names (Reaction)
        test_object = cr.create_object(
            identifier="AMONITRO-RXN",
            directory=dir_data,
            compartment="c",
            database="META",
            replacement=test_replacement,
        )
        self.assertIsInstance(obj=test_object, cls=cobra_core.Reaction)
        self.assertCountEqual(
            first=getattr(test_object, "compartments"), second={"c"}
        )
        self.assertEqual(
            first=getattr(test_object, "id"), second="monooxygenase_c"
        )
        self.assertIn(
            member="Donor_new_name_c",
            container=[meta.id for meta in getattr(test_object, "metabolites")],
        )

    def test_add_metabolites(self):
        # CASE: From path
        test_model = cobra_core.Model(NAME)
        cr.add_metabolites(
            model=test_model,
            obj=dir_input.joinpath("metabolites_01_normal.txt"),
            directory=dir_data,
            database="META",
        )
        # CASE: path as str
        cr.add_metabolites(
            model=test_model,
            obj=str(dir_input.joinpath("metabolites_01_normal.txt")),
            directory=dir_data,
            database="META",
        )
        test_names = ["HOMOMETHIONINE_c", "MALTOSE_b"]
        self.assertListEqual(
            list1=[member.id for member in test_model.metabolites],
            list2=test_names,
        )
        # CASE: From string
        test_model = cobra_core.Model(NAME)
        test_string = ["HOMOMETHIONINE, c"]
        cr.add_metabolites(
            model=test_model,
            obj=test_string,
            directory=dir_data,
            database="META",
        )
        self.assertIn(
            member="HOMOMETHIONINE_c",
            container=[member.id for member in test_model.metabolites],
        )
        # CASE: From List, Custom
        test_model = cobra_core.Model(NAME)
        test_string = ["Custom_c, Custom metabolite, c, H20, 0 "]
        cr.add_metabolites(
            model=test_model,
            obj=test_string,
            directory=dir_data,
            database=None,
        )
        self.assertIn(
            member="Custom_c",
            container=[member.id for member in test_model.metabolites],
        )
        # CASE: multiple entries
        test_model = cobra_core.Model(NAME)
        test_list = ["HOMOMETHIONINE, c", "MALTOSE, c"]
        cr.add_metabolites(
            model=test_model,
            obj=test_list,
            directory=dir_data,
            database="META",
        )
        test_names = ["HOMOMETHIONINE_c", "MALTOSE_c"]
        self.assertListEqual(
            list1=[member.id for member in test_model.metabolites],
            list2=test_names,
        )
        # CASE: Multiple metabolites objects
        test_model = cobra_core.Model(NAME)
        test_list = [
            textbook_kegg.metabolites.get_by_id(item)
            for item in ("C00001_c", "C00002_c", "C00003_c")
        ]
        cr.add_metabolites(model=test_model, obj=test_list)
        for item in ("C00001_c", "C00002_c", "C00003_c"):
            self.assertIn(
                member=item,
                container=[member.id for member in test_model.metabolites],
            )
        # CASE: BIGG database
        test_model = cobra_core.Model(NAME)
        cr.add_metabolites(
            model=test_model,
            obj=["acald, c"],
            directory=dir_data,
            database="BIGG",
            model_id="universal",
        )
        self.assertIn(
            member="acald_c",
            container=[member.id for member in test_model.metabolites],
        )
        # CASE: Using replacement dictionary
        test_model = cobra_core.Model(NAME)
        test_list = ["HOMOMETHIONINE, c", "MALTOSE, c"]
        cr.add_metabolites(
            model=test_model,
            obj=test_list,
            directory=dir_data,
            database="META",
            replacement={"HOMOMETHIONINE": "WATER"},
        )
        self.assertIn(
            member="WATER_c",
            container=[member.id for member in test_model.metabolites],
        )
        self.assertNotIn(
            member="HOMOMETHIONINE_c",
            container=[member.id for member in test_model.metabolites],
        )
        self.assertIn(
            member="MALTOSE_c",
            container=[member.id for member in test_model.metabolites],
        )

    def test_add_reactions(self):
        # CASE: From str
        test_model = cobra_core.Model(NAME)
        cr.add_reactions(
            model=test_model,
            obj=str(dir_input.joinpath("reactions_normal.txt")),
            directory=dir_data,
            database="META",
        )
        for reaction in ("GLC_cb", "RXN_14462_p"):
            self.assertIn(
                member=reaction,
                container=[reaction.id for reaction in test_model.reactions],
            )
        # CASE: From Path
        test_model = cobra_core.Model(NAME)
        cr.add_reactions(
            model=test_model,
            obj=dir_input.joinpath("reactions_normal.txt"),
            directory=dir_data,
            database="META",
        )
        for reaction in ("GLC_cb", "RXN_14462_p"):
            self.assertIn(
                member=reaction,
                container=[reaction.id for reaction in test_model.reactions],
            )
        # CASE: Custom reaction
        test_model = cobra_core.Model(NAME)
        cr.add_reactions(
            model=test_model,
            obj=["GLC_cb, Glucose Transport|GLC_c <-> GLC_b"],
            directory=dir_data,
            database="META",
        )
        self.assertIn(
            member="GLC_cb",
            container=[reaction.id for reaction in test_model.reactions],
        )

        # CASE: From string, custom metabolites
        test_model = cobra_core.Model(NAME)
        cr.add_reactions(
            model=test_model,
            obj=["Custom_cb, Custom reaction|Custom_c <=> Custom_b"],
            directory=dir_data,
            database=None,
        )
        self.assertIn(
            member="Custom_cb",
            container=[reaction.id for reaction in test_model.reactions],
        )
        for metabolite in ["Custom_c", "Custom_b"]:
            self.assertIn(
                member=metabolite,
                container=[meta.id for meta in test_model.metabolites],
            )
        # CASE: From string, KEGG
        test_model = cobra_core.Model(NAME)
        cr.add_reactions(
            model=test_model,
            obj=["R02736, c"],
            directory=dir_data,
            database="KEGG",
            genome="hsa",
        )
        self.assertIn(
            member="R02736_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE: From List of strings
        test_model = cobra_core.Model(NAME)
        test_list = [
            "GLC_cb, Glucose Transport|GLC_c <-> GLC_b",
            "RXN-14462, c",
        ]
        cr.add_reactions(
            model=test_model,
            obj=test_list,
            directory=dir_data,
            database="GCF_000020025",
        )
        for reaction in ("GLC_cb", "RXN_14462_c"):
            self.assertIn(
                member=reaction,
                container=[reaction.id for reaction in test_model.reactions],
            )
        self.assertIn(
            member="NPUN_RS12370",
            container=[gene.id for gene in test_model.genes],
        )

        # CASE: multiple reactions
        test_model = cobra_core.Model(NAME)
        test_list = [
            textbook_kegg.reactions.get_by_id(reaction)
            for reaction in ("ACALDt", "ATPS4r", "ACt2r")
        ]
        cr.add_reactions(model=test_model, obj=test_list, directory=dir_data)
        for reaction in ("ACALDt", "ATPS4r", "ACt2r"):
            self.assertIn(
                member=reaction,
                container=[reaction.id for reaction in test_model.reactions],
            )

        # CASE: BIGG
        test_model = cobra_core.Model(NAME)
        cr.add_reactions(
            model=test_model,
            obj=["ACALDt, c"],
            directory=dir_data,
            database="BIGG",
            model_id="universal",
        )
        self.assertIn(
            member="ACALDt_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE: Duplicate element
        test_model = textbook_kegg
        cr.add_reactions(
            model=test_model,
            obj=["R08549, c"],
            directory=dir_data,
            database="KEGG",
        )
        cr.add_reactions(
            model=test_model,
            obj=["AKGDH, c"],
            directory=dir_data,
            database="BIGG",
            model_id="universal",
        )
        self.assertNotIn(
            member="AKGDH_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE: replacement dicts
        test_model = cobra_core.Model(NAME)
        test_list = [
            "GLC_cb, Glucose Transport|Not_GLC_c <-> GLC_b",
            "RXN-14462, c",
            "ACETALD-DEHYDROG-RXN ,c",
            "AMONITRO-RXN, c",
        ]
        cr.add_reactions(
            model=test_model,
            obj=test_list,
            directory=dir_data,
            database="META",
            replacement={
                "Not_GLC": "GLC",
                "ACETALD-DEHYDROG-RXN": "ACETALD-DEHYDROG-RXN-NEW",
                "Donor-H2": "Donor_Water",
            },
        )
        self.assertTrue(
            test_model.reactions.has_id("GLC_cb"),
        )
        self.assertTrue(
            test_model.reactions.has_id("ACETALD_DEHYDROG_RXN_NEW_c"),
        )
        self.assertFalse(
            test_model.reactions.has_id("Not_GLC_c"),
        )
        for metabolite in ("GLC_b", "GLC_c", "Donor_Water_c"):
            self.assertTrue(
                test_model.metabolites.has_id(metabolite),
                msg=f"{metabolite} not found",
            )

        test_model = cobra_core.Model(NAME)
        test_list = [
            "GLC_cb, Glucose Transport| GLC_c <-> GLC_b",
            "RXN_17742_c, RXN_17742_c |"
            + "1 Oxidized-ferredoxins_c <-> 1 Reduced-ferredoxins_c ",
        ]
        cr.add_reactions(model=test_model, obj=test_list, directory=dir_data)

        for meta in ("GLC_cb", "RXN_17742_c"):
            self.assertTrue(test_model.reactions.has_id(meta))

        # CASE PMN:ARA
        test_model = cobra_core.Model(NAME)
        reaction_str = [
            "PYR_MAL_pc, Pyruvate/Malate transporter |"
            "PYRUVATE_p + MAL_c <-> PYRUVATE_c + MAL_p"
        ]

        cr.add_reactions(test_model, reaction_str, dir_data, "PMN:ARA")

        self.assertTrue(test_model.reactions.has_id("PYR_MAL_pc"))


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
