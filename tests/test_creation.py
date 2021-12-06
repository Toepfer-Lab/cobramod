#!/usr/bin/env python3
"""Unittest for module creation

In this test, multiple single functions for the creation of reactions and
metabolites are tested. Important is the independency of the functions, unless
they are part of the API such as add_metabolites. Therefore the modules is
separated in:

- SimpleFunctions: Creation of objects
- ComplexFunctions: Functions, that uses multiple simple functions.
"""
from logging import DEBUG
from pathlib import Path
from unittest import main, TestCase

from cobra import Metabolite, Reaction, Model, __version__
from requests import HTTPError

from cobramod.core import creation as cr
from cobramod.core.retrieval import get_data
from cobramod.debug import debug_log
from cobramod.error import WrongSyntax, NoDelimiter
from cobramod.parsing import BaseParser
from cobramod.test import textbook_kegg
from cobramod import __version__ as cobramod_version

# Debug must be set in level DEBUG for the test
debug_log.setLevel(DEBUG)
# Setting directory for datA
dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")
print(f"CobraMod version: {cobramod_version}")
print(f"COBRApy version: {__version__}")


class SimpleFunctions(TestCase):
    """
    Test for simple test such as creating metabolites from string or from
    files.
    """

    # TODO: use replacement dictionaries !!

    def test__metabolite_from_string(self):
        # CASE 1: Correct input
        test_string = "MALTOSE_b, MALTOSE[b], b, C12H22O11, 0"
        test_metabolite = cr._metabolite_from_string(line_string=test_string)
        # Checking that new meta is in Model
        self.assertIsInstance(obj=test_metabolite, cls=Metabolite)
        self.assertEqual(first=test_metabolite.id, second="MALTOSE_b")
        self.assertEqual(first=test_metabolite.name, second="MALTOSE[b]")
        self.assertEqual(first=test_metabolite.charge, second=0)
        self.assertDictEqual(
            d1=test_metabolite.elements, d2={"C": 12, "H": 22, "O": 11}
        )

    def test__get_metabolite(self):
        # CASE 1: regular META
        test_dict = get_data(
            directory=dir_data,
            identifier="HOMOMETHIONINE",
            database="META",
            debug_level=10,
        )
        test_metabolite = cr._get_metabolite(
            metabolite_dict=test_dict, compartment="c"
        )
        self.assertIsInstance(obj=test_metabolite, cls=Metabolite)
        self.assertEqual(first=test_metabolite.id, second="HOMOMETHIONINE_c")
        self.assertEqual(first=test_metabolite.name, second="L-homomethionine")
        self.assertEqual(first=test_metabolite.formula, second="C6H13N1O2S1")
        # CASE 2a: Find translation
        test_dict = get_data(
            directory=dir_data,
            identifier="WATER",
            database="META",
            debug_level=10,
        )
        test_metabolite = cr._get_metabolite(
            metabolite_dict=test_dict, compartment="c", model=textbook_kegg
        )
        self.assertEqual(first=test_metabolite.id, second="C00001_c")
        # CASE 2b: Find translation (BIGG)
        test_dict = get_data(
            directory=dir_data,
            identifier="h2o",
            database="BIGG",
            debug_level=10,
            model_id="universal",
        )
        test_metabolite = cr._get_metabolite(
            metabolite_dict=test_dict, compartment="c", model=textbook_kegg
        )
        self.assertEqual(first=test_metabolite.id, second="C00001_c")
        # TODO: add extra cases, ARA, KEGG

    def test__convert_string_metabolite(self):
        # CASE 1: retrieval from META
        test_metabolite = cr._convert_string_metabolite(
            line="HOMOMETHIONINE, c",
            model=Model(0),
            directory=dir_data,
            database="META",
        )
        self.assertIsInstance(obj=test_metabolite, cls=Metabolite)
        self.assertEqual(first=test_metabolite.id, second="HOMOMETHIONINE_c")
        self.assertEqual(first=test_metabolite.name, second="L-homomethionine")
        self.assertEqual(first=test_metabolite.formula, second="C6H13N1O2S1")
        self.assertEqual(first=test_metabolite.compartment, second="c")
        # CASE 2: custom metabolite
        test_metabolite = cr._convert_string_metabolite(
            line="MALTOSE_b, MALTOSE[b], b, C12H22O11, 0",
            model=Model(0),
            directory=dir_data,
        )
        self.assertEqual(first=test_metabolite.id, second="MALTOSE_b")
        self.assertEqual(first=test_metabolite.name, second="MALTOSE[b]")
        self.assertEqual(first=test_metabolite.charge, second=0)
        self.assertDictEqual(
            d1=test_metabolite.elements, d2={"C": 12, "H": 22, "O": 11}
        )

    def test__get_file_metabolites(self):
        test_model = Model(0)
        # CASE 0: Testing if file is not found
        self.assertRaises(
            FileNotFoundError,
            cr._get_file_metabolites,
            test_model,
            Path.cwd().joinpath("no_file"),
        )
        # CASE 1: Metabolite is not found (or misspelled)
        self.assertRaises(
            HTTPError,
            cr._get_file_metabolites,
            model=test_model,
            filename=dir_input.joinpath("metabolites_02_misspelled.txt"),
            # Directory to save / check for xml files
            directory=dir_data,
            database="META",
        )
        # CASE 2: Bad format (e. g. charge is missing).
        self.assertRaises(
            WrongSyntax,
            cr._get_file_metabolites,
            model=test_model,
            filename=dir_input.joinpath("metabolites_03_badFormat.txt"),
            # Directory to save / check for xml files
            directory=dir_data,
            database="META",
        )
        # CASE 3: Normal input
        test_list = cr._get_file_metabolites(
            model=test_model,
            # File with Metabolites to add
            filename=dir_input.joinpath("metabolites_01_normal.txt"),
            # Directory to save / check for xml files
            directory=dir_data,
            database="META",
        )
        # Length should increase
        self.assertEqual(len(test_list), 2)
        # Both Names in metabolites
        test_names = ["HOMOMETHIONINE_c", "MALTOSE_b"]
        self.assertListEqual(
            list1=[member.id for member in test_list], list2=test_names
        )

    def test__get_reaction(self):
        # CASE 1: Regular Biocyc reaction
        test_data = get_data(
            directory=dir_data,
            identifier="OXALODECARB-RXN",
            database="ARA",
            debug_level=10,
        )
        test_reaction = cr._get_reaction(
            data_dict=test_data,
            compartment="c",
            directory=dir_data,
            database="ARA",
            replacement={},
            model=Model(0),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIsInstance(obj=test_reaction, cls=Reaction)
        self.assertTupleEqual(tuple1=test_reaction.bounds, tuple2=(0, 1000))
        test_list = [
            "PYRUVATE_c",
            "CARBON_DIOXIDE_c",
            "OXALACETIC_ACID_c",
            "PROTON_c",
        ]
        self.assertListEqual(
            list1=sorted([meta.id for meta in test_reaction.metabolites]),
            list2=sorted(test_list),
        )
        # CASE 2: Regular KEGG reaction in compartment "p". Irreversible
        test_data = get_data(
            directory=dir_data,
            identifier="R02736",
            database="KEGG",
            debug_level=10,
            genome="hsa",
        )
        test_reaction = cr._get_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="KEGG",
            replacement={},
            model=Model(0),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIsInstance(obj=test_reaction, cls=Reaction)
        self.assertTupleEqual(tuple1=test_reaction.bounds, tuple2=(0, 1000))
        test_list = [
            "C01172_p",
            "C00006_p",
            "C01236_p",
            "C00005_p",
            "C00080_p",
        ]
        self.assertListEqual(
            list1=sorted([meta.id for meta in test_reaction.metabolites]),
            list2=sorted(test_list),
        )
        # CASE 3: Reactions with different coefficients
        test_data = get_data(
            directory=dir_data,
            identifier="R00114",
            database="KEGG",
            debug_level=10,
            genome="eco",
        )
        test_reaction = cr._get_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="KEGG",
            replacement={},
            model=Model(0),
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertIsInstance(obj=test_reaction, cls=Reaction)
        self.assertTupleEqual(
            tuple1=test_reaction.bounds, tuple2=(-1000, 1000)
        )
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
        # CASE 4: transport reaction, simple
        test_data = get_data(
            directory=dir_data,
            identifier="TRANS-RXN0-574",
            database="META",
            debug_level=10,
        )
        test_reaction = cr._get_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="META",
            replacement={},
            model=Model(0),
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
        # CASE 5: transport reaction, complex
        test_data = get_data(
            directory=dir_data,
            identifier="ABC-56-RXN",
            database="META",
            debug_level=10,
        )
        test_reaction = cr._get_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="META",
            replacement={},
            model=Model(0),
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
        # CASE 6a: Retrieve reaction with translated equivalents
        test_data = get_data(
            directory=dir_data,
            identifier="ACALD",
            database="BIGG",
            debug_level=10,
            model_id="e_coli_core",
        )
        test_reaction = cr._get_reaction(
            data_dict=test_data,
            compartment="c",
            directory=dir_data,
            database="BIGG",
            replacement={},
            model=textbook_kegg,
            model_id="universal",
            show_imbalance=True,
            stop_imbalance=False,
        )
        self.assertEqual(first=test_reaction.id, second="R00228_c")
        # CASE 6b: Retrieve reaction with translated equivalents. Check for
        # metabolites.
        test_data = get_data(
            directory=dir_data,
            identifier="ADENODEAMIN-RXN",
            database="META",
            debug_level=10,
        )
        # FIXME: find a solution for unformatted identifiers
        test_reaction = cr._get_reaction(
            data_dict=test_data,
            compartment="c",
            directory=dir_data,
            database="META",
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

    def test__convert_string_reaction(self):
        # CASE 1: using delimiter, compartment is cytosol
        test_model = Model(0)
        test_line = (
            "RXN_17742_c, RXN_17742_c |"
            "1 Oxidized-ferredoxins_c <-> 1 Reduced-ferredoxins_c "
        )
        test_reaction = cr._convert_string_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
        )
        self.assertEqual(first="RXN_17742_c", second=test_reaction.id)
        # CASE 2: No delimiter
        test_model = Model(0)
        test_line = "RXN-14462, c"
        test_reaction = cr._convert_string_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
        )
        self.assertEqual(first="RXN_14462_c", second=test_reaction.id)
        # CASE 2: No delimiter, compartment p
        test_model = Model(0)
        test_line = "RXN-14462, p"
        test_reaction = cr._convert_string_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="GCF_000020025",
            stop_imbalance=False,
            show_imbalance=True,
        )
        self.assertEqual(first="RXN_14462_p", second=test_reaction.id)
        # CASE 3: 100% Custom metabolite
        test_model = Model(0)
        test_line = "CUSTOM_rxn1_p, Custom_reaction | Meta_A_p --> Meta_B_p"
        test_reaction = cr._convert_string_reaction(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
        )
        self.assertEqual(first="CUSTOM_rxn1_p", second=test_reaction.id)
        for metabolite in ["Meta_A_p", "Meta_B_p"]:
            self.assertIn(
                member=metabolite,
                container=[meta.id for meta in test_reaction.metabolites],
            )

    def test__obtain_reaction(self):
        # CASE 1: Regular META reaction
        test_model = Model(0)
        test_model.compartments = {"e": "extracellular", "p": "plastid"}
        test_reaction = cr._obtain_reaction(
            model=test_model,
            directory=dir_data,
            identifier="OXALODECARB-RXN",
            database="VCHO",
            compartment="p",
            replacement={},
            show_imbalance=True,
            stop_imbalance=False,
            genome=None,
            model_id=None,
        )
        self.assertEqual(first="OXALODECARB_RXN_p", second=test_reaction.id)
        self.assertCountEqual(
            first=[gene.id for gene in test_reaction.genes],
            second=["VC0550", "VC0551", "VC0792"],
        )
        # CASE 2: check for equivalent. (Similar to CASE 6b in _get_reaction)
        test_model = textbook_kegg.copy()
        test_reaction = cr._obtain_reaction(
            model=test_model,
            directory=dir_data,
            compartment="c",
            database="META",
            replacement={},
            identifier="ADENODEAMIN-RXN",
            show_imbalance=True,
            stop_imbalance=False,
            genome=None,
            model_id=None,
        )
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
        # TODO: test more databases

    def test__reaction_information(self):
        # CASE 0: Missing compartment
        self.assertRaises(
            WrongSyntax, cr._reaction_information, string=["2 GLC_c <=> GLC"]
        )
        # CASE 1: Regular usage
        test_tuple, test_dict = cr._reaction_information(
            string="GLC_c --> GLC_b"
        )
        self.assertEqual(first=test_tuple, second=(0, 1000))
        self.assertDictEqual(d1=test_dict, d2={"GLC_c": -1.0, "GLC_b": 1.0})
        # CASE 2: Regular usage. One-sided
        test_tuple, test_dict = cr._reaction_information(string="GLC_c -->")
        self.assertEqual(first=test_tuple, second=(0, 1000))
        self.assertDictEqual(d1=test_dict, d2={"GLC_c": -1.0})

    def test__reaction_from_string(self):
        # CASE 1: wrong format, no delimiter
        self.assertRaises(
            NoDelimiter,
            cr._reaction_from_string,
            line_string="GLC_cb, GLC_cb GLC_c <-> GLC_b",
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
        )
        # CASE 2: Normal, ID and name differ
        test_reaction = cr._reaction_from_string(
            line_string="GLC_cb, Glucose Transport| 2 GLC_c <=> GLC_b",
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=False,
        )
        # Checking if ID, name and instance are correct.
        self.assertEqual(first=test_reaction.id, second="GLC_cb")
        self.assertEqual(first=test_reaction.name, second="Glucose Transport")
        self.assertIsInstance(obj=test_reaction, cls=Reaction)
        test_list = ["GLC_c", "GLC_b"]
        # All metabolites should be there
        for metabolite in test_list:
            self.assertIn(
                member=metabolite,
                container=[meta.id for meta in test_reaction.metabolites],
            )
        self.assertEqual(
            first=-2, second=test_reaction.get_coefficient("GLC_c")
        )
        self.assertEqual(
            first=1, second=test_reaction.get_coefficient("GLC_b")
        )
        # CASE 3: Larger reaction
        test_reaction = cr._reaction_from_string(
            line_string="GLC_cb, Glucose Transport|"
            + "4 GLT_c + 2 GLC_c <=> GLC_b + 2 GLT_b",
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=False,
        )
        self.assertEqual(
            first=-4, second=test_reaction.get_coefficient("GLT_c")
        )
        self.assertEqual(
            first=2, second=test_reaction.get_coefficient("GLT_b")
        )

    def test__get_file_reactions(self):
        test_model = Model(0)
        test_list = cr._get_file_reactions(
            model=test_model,
            filename=dir_input.joinpath("reactions_normal.txt"),
            directory=dir_data,
            database="META",
            stop_imbalance=False,
            show_imbalance=True,
        )
        test_names = ["GLC_cb", "RXN_14462_p"]
        self.assertListEqual(
            list1=test_names, list2=[reaction.id for reaction in test_list]
        )


class ComplexFunctions(TestCase):
    """
    Test of functions for the API
    """

    @classmethod
    def setUpClass(cls):
        BaseParser.ignore_db_versions = True

    def test_create_object(self):
        # CASE 1a: metabolite from MetaCyc
        test_object = cr.create_object(
            identifier="GLC",
            directory=dir_data,
            compartment="c",
            database="META",
        )
        self.assertIsInstance(obj=test_object, cls=Metabolite)
        self.assertEqual(first=test_object.compartment, second="c")
        # CASE 1b: metabolite from KEGG
        test_object = cr.create_object(
            identifier="C00026",
            directory=dir_data,
            compartment="p",
            database="KEGG",
        )
        self.assertIsInstance(obj=test_object, cls=Metabolite)
        self.assertEqual(first=test_object.compartment, second="p")
        # CASE 2a: reaction from MetaCyc
        test_object = cr.create_object(
            identifier="2.6.1.58-RXN",
            directory=dir_data,
            compartment="p",
            database="META",
        )
        self.assertIsInstance(obj=test_object, cls=Reaction)
        self.assertIn(member="p", container=test_object.compartments)
        # CASE 2b: reaction from KEGG
        test_object = cr.create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
            genome="hsa",
        )
        self.assertIsInstance(obj=test_object, cls=Reaction)
        self.assertIn(member="c", container=test_object.compartments)
        self.assertCountEqual(
            first=[gene.id for gene in test_object.genes],
            second=["9563", "2539"],
        )
        # CASE 2c: Reaction of Biocyc which could be a pathway as well
        test_object = cr.create_object(
            identifier="AMONITRO-RXN",
            directory=dir_data,
            compartment="c",
            database="META",
            show_imbalance=False,
        )
        self.assertIsInstance(obj=test_object, cls=Reaction)
        # CASE 3a: pathway from MetaCyc
        test_object = cr.create_object(
            identifier="PWY-1187",
            directory=dir_data,
            compartment="c",
            database="META",
        )
        self.assertIsInstance(obj=test_object, cls=dict)
        self.assertEqual(first=test_object["ENTRY"], second="PWY-1187")
        # CASE 3a: pathway from MetaCyc
        test_object = cr.create_object(
            identifier="M00001",
            directory=dir_data,
            compartment="c",
            database="KEGG",
        )
        self.assertIsInstance(obj=test_object, cls=dict)
        self.assertEqual(first=test_object["ENTRY"], second="M00001")

    def test_add_metabolites(self):
        # CASE 0: Missing Arguments
        self.assertRaises(
            ValueError,
            cr.add_metabolites,
            model=Model(0),
            obj=dir_input.joinpath("metabolites_01_normal.txt"),
        )
        # CASE 1: From path
        test_model = Model(0)
        cr.add_metabolites(
            model=test_model,
            obj=dir_input.joinpath("metabolites_01_normal.txt"),
            directory=dir_data,
            database="META",
        )
        test_names = ["HOMOMETHIONINE_c", "MALTOSE_b"]
        self.assertListEqual(
            list1=[member.id for member in test_model.metabolites],
            list2=test_names,
        )
        # CASE 2: From string
        test_model = Model(0)
        test_string = "HOMOMETHIONINE, c"
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
        # CASE 2: From string, Custom
        test_model = Model(0)
        test_string = "Custom_c, Custom metabolites, c, H20, 0 "
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
        # CASE 3: From List of strings
        test_model = Model(0)
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
        # CASE 4: In case of single metabolite
        test_model = Model(0)
        test_metabolite = textbook_kegg.metabolites.get_by_id("C00001_c")
        cr.add_metabolites(model=test_model, obj=test_metabolite)
        self.assertIn(
            member="C00001_c",
            container=[member.id for member in test_model.metabolites],
        )
        # CASE 5: In case of multiple metabolites
        test_model = Model(0)
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
        # CASE 6: CASE BIGG
        test_model = Model(0)
        cr.add_metabolites(
            model=test_model,
            obj="acald, c",
            directory=dir_data,
            database="BIGG",
            model_id="universal",
        )
        self.assertIn(
            member="acald_c",
            container=[member.id for member in test_model.metabolites],
        )

    def test_add_reactions(self):
        # CASE 0: Missing arguments.
        self.assertRaises(
            ValueError,
            cr.add_reactions,
            model=Model(0),
            obj=dir_input.joinpath("reactions_normal.txt"),
        )
        # CASE 1: From Path
        test_model = Model(0)
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
        # CASE 2a: From string
        test_model = Model(0)
        cr.add_reactions(
            model=test_model,
            obj="GLC_cb, Glucose Transport|GLC_c <-> GLC_b",
            directory=dir_data,
            database="META",
        )
        self.assertIn(
            member="GLC_cb",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE 2b: From string, custom metabolites
        test_model = Model(0)
        cr.add_reactions(
            model=test_model,
            obj="Custom_cb, Custom reaction|Custom_c <=> Custom_b",
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
        # CASE 2c: From string, KEGG
        test_model = Model(0)
        cr.add_reactions(
            model=test_model,
            obj="R02736, c",
            directory=dir_data,
            database="KEGG",
            genome="hsa",
        )
        self.assertIn(
            member="R02736_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE 3: From List of strings
        test_model = Model(0)
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
        # CASE 4: In case of single reaction
        test_model = Model(0)
        test_reaction = textbook_kegg.reactions.get_by_id("ACALDt")
        cr.add_reactions(model=test_model, obj=test_reaction)
        self.assertIn(
            member="ACALDt",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE 5: In case of multiple reactions
        test_model = Model(0)
        test_list = [
            textbook_kegg.reactions.get_by_id(reaction)
            for reaction in ("ACALDt", "ATPS4r", "ACt2r")
        ]
        cr.add_reactions(model=test_model, obj=test_list)
        for reaction in ("ACALDt", "ATPS4r", "ACt2r"):
            self.assertIn(
                member=reaction,
                container=[reaction.id for reaction in test_model.reactions],
            )
        # CASE 5: BIGG
        test_model = Model(0)
        cr.add_reactions(
            model=test_model,
            obj="ACALDt, c",
            directory=dir_data,
            database="BIGG",
            model_id="universal",
        )
        self.assertIn(
            member="ACALDt_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE 6: Duplicate element
        test_model = textbook_kegg
        cr.add_reactions(
            model=test_model,
            obj="R08549, c",
            directory=dir_data,
            database="KEGG",
        )
        cr.add_reactions(
            model=test_model,
            obj="AKGDH, c",
            directory=dir_data,
            database="BIGG",
            model_id="universal",
        )
        self.assertNotIn(
            member="AKGDH_c",
            container=[reaction.id for reaction in test_model.reactions],
        )


if __name__ == "__main__":
    main(verbosity=2)
