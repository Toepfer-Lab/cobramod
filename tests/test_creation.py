#!/usr/bin/env python3
from pathlib import Path
from cobramod import creation as cr
from cobramod.mod_parser import get_data
from cobramod.debug import debug_log
from logging import DEBUG
import unittest
import cobra as cb

# configuration
debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")

if not dir_data.exists():
    dir_data.mkdir(parents=True)


class ModulTesting(unittest.TestCase):
    # TODO: use replacement dictionaries !!

    def test__create_meta_from_string(self):
        # CASE 1: Correct input
        testInput = "MALTOSE_b, MALTOSE[b], b, C12H22O11, 0"
        testMeta = cr._create_meta_from_string(line_string=testInput)
        # Checking that new meta is in Model
        self.assertIsInstance(testMeta, cb.Metabolite)
        # TODO: verify attributes

    def test_build_metabolite(self):
        # CASE 1: regular META
        test_dict = get_data(
            directory=dir_data,
            identifier="HOMOMETHIONINE",
            database="META",
            debug_level=10,
        )
        test_metabolite = cr.build_metabolite(
            metabolite_dict=test_dict, compartment="c"
        )
        self.assertIsInstance(obj=test_metabolite, cls=cb.Metabolite)
        self.assertEqual(first=test_metabolite.id, second="HOMOMETHIONINE_c")
        self.assertEqual(first=test_metabolite.name, second="L-homomethionine")
        self.assertEqual(first=test_metabolite.formula, second="C6H13N1O2S1")
        # TODO: add extra cases, ARA, KEGG

    def test_meta_string_to_model(self):
        # CASE 1: retrieval from META
        test_metabolite = cr.meta_string_to_model(
            line="HOMOMETHIONINE, c",
            model=cb.Model(0),
            directory=dir_data,
            database="META",
        )
        self.assertIsInstance(obj=test_metabolite, cls=cb.Metabolite)
        self.assertEqual(first=test_metabolite.id, second="HOMOMETHIONINE_c")
        self.assertEqual(first=test_metabolite.name, second="L-homomethionine")
        self.assertEqual(first=test_metabolite.formula, second="C6H13N1O2S1")
        self.assertEqual(first=test_metabolite.compartment, second="c")
        # TODO: extra cases
        # CASE 2: custom metabolite

    def test_add_meta_from_file(self):
        test_model = cb.Model(0)
        # Testing if model is not cobra.Model
        self.assertRaises(
            TypeError,
            cr.add_meta_from_file,
            "notModel",
            Path.cwd().joinpath("nofile"),
        )
        # Testing if file is not found
        self.assertRaises(
            FileNotFoundError,
            cr.add_meta_from_file,
            test_model,
            Path.cwd().joinpath("nofile"),
        )
        # CASE 1: Metabolite is not found (or misspelled)
        self.assertRaises(
            Warning,
            cr.add_meta_from_file,
            model=test_model,
            filename=dir_input.joinpath("metaToAdd_02_misspelled.txt"),
            # Directory to save / check for xml files
            directory=dir_data,
            database="META",
        )
        # CASE 2: Bad format (e. g. charge is missing).
        self.assertRaises(
            IndexError,
            cr.add_meta_from_file,
            model=test_model,
            filename=dir_input.joinpath("metaToAdd_03_badFormat.txt"),
            # Directory to save / check for xml files
            directory=dir_data,
            database="META",
        )
        # CASE 3: Normal input
        cr.add_meta_from_file(
            model=test_model,
            # File with Metabolites to add
            filename=dir_input.joinpath("metaToAdd_01_normal.txt"),
            # Directory to save / check for xml files
            directory=dir_data,
            database="META",
        )
        # Length should increase
        self.assertEqual(len(test_model.metabolites), 2)
        # Both Names in metabolites
        test_list = ["HOMOMETHIONINE_c", "MALTOSE_b"]
        self.assertTrue(
            all(
                [
                    x in [meta.id for meta in test_model.metabolites]
                    for x in test_list
                ]
            )
        )

    def test__build_reaction(self):
        # CASE 1: Regular Biocyc reaction
        test_data = get_data(
            directory=dir_data,
            identifier="OXALODECARB-RXN",
            database="META",
            debug_level=10,
        )
        test_reaction = cr._build_reaction(
            data_dict=test_data,
            compartment="c",
            directory=dir_data,
            database="META",
            replacement_dict={},
        )
        self.assertIsInstance(obj=test_reaction, cls=cb.Reaction)
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
        # CASE 2: Regular KEGG reaction in compartement "p"
        test_data = get_data(
            directory=dir_data,
            identifier="R02736",
            database="KEGG",
            debug_level=10,
        )
        test_reaction = cr._build_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="KEGG",
            replacement_dict={},
        )
        self.assertIsInstance(obj=test_reaction, cls=cb.Reaction)
        self.assertTupleEqual(
            tuple1=test_reaction.bounds, tuple2=(-1000, 1000)
        )
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
        )
        test_reaction = cr._build_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="KEGG",
            replacement_dict={},
        )
        self.assertIsInstance(obj=test_reaction, cls=cb.Reaction)
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
        test_reaction = cr._build_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="META",
            replacement_dict={},
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
        test_reaction = cr._build_reaction(
            data_dict=test_data,
            compartment="p",
            directory=dir_data,
            database="META",
            replacement_dict={},
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

    def test__add_reaction_line_to_model(self):
        # CASE 1: using delimiter, compartment is cytosol
        test_model = cb.Model(0)
        test_line = (
            "RXN_17742_c, RXN_17742_c |"
            "Oxidized-ferredoxins_c:-1, Reduced-ferredoxins_c: 1"
        )
        cr._add_reaction_line_to_model(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
        )
        self.assertTrue(
            "RXN_17742_c" in [reaction.id for reaction in test_model.reactions]
        )
        # CASE 2: No delimiter
        test_model = cb.Model(0)
        test_line = "RXN-14462, c"
        cr._add_reaction_line_to_model(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
        )
        self.assertTrue(
            "RXN_14462_c" in [reaction.id for reaction in test_model.reactions]
        )
        # CASE 2: No delimiter, compartment p
        test_model = cb.Model(0)
        test_line = "RXN-14462, p"
        cr._add_reaction_line_to_model(
            line=test_line,
            model=test_model,
            directory=dir_data,
            database="META",
        )
        self.assertTrue(
            "RXN_14462_p" in [reaction.id for reaction in test_model.reactions]
        )

    def test_add_reaction(self):
        # CASE 1: Regular META reaction
        test_model = cb.Model(0)
        test_model.compartments = {"e": "extracellular", "p": "plastid"}
        cr.add_reaction(
            model=test_model,
            directory=dir_data,
            identifier="OXALODECARB-RXN",
            database="META",
            compartment="p",
            replacement_dict={},
        )
        self.assertTrue(
            "OXALODECARB_RXN_p" in [rxn.id for rxn in test_model.reactions]
        )

    def test__build_dict_for_metabolites(self):
        # CASE 0a: TypeError
        self.assertRaises(TypeError, cr._build_dict_for_metabolites, str())
        # CASE 0b: Coefficient missing
        self.assertRaises(
            ValueError, cr._build_dict_for_metabolites, string_list=[" GLC_c:"]
        )
        self.assertDictEqual(
            {"GLC_b": 1.0, "GLC_c": -1.0},
            cr._build_dict_for_metabolites(
                string_list=[" GLC_c:-1", "GLC_b: 1"]
            ),
        )

    def test_create_custom_reaction(self):
        # CASE 0: wrong format, no delimiter
        self.assertRaises(
            IndexError,
            cr.create_custom_reaction,
            line_string="GLC_cb, GLC_cb GLC_c:-1, GLC_b:1",
            directory=dir_data,
            database="META",
        )
        # CASE 1: No ID detected
        self.assertRaises(
            Warning,
            cr.create_custom_reaction,
            line_string=" |GLC_c:-1, GLC_b:1",
            directory=dir_data,
            database="META",
        )
        # CASE 2: Normal, ID and name differ
        ReactionTest = cr.create_custom_reaction(
            line_string="GLC_cb, Glucose Transport|GLC_c:-1, GLC_b:1",
            directory=dir_data,
            database="META",
        )
        # Checking if ID, name and instance are correct.
        self.assertTrue(
            all(
                [
                    ReactionTest.id == "GLC_cb",
                    ReactionTest.name == "Glucose Transport",
                    isinstance(ReactionTest, cb.Reaction),
                ]
            )
        )
        test_list = ["GLC_c", "GLC_b"]
        # names should be there
        self.assertTrue(
            all(
                [
                    x in [meta.id for meta in ReactionTest.metabolites]
                    for x in test_list
                ]
            )
        )
        self.assertEqual(-1, ReactionTest.get_coefficient("GLC_c"))
        self.assertEqual(1, ReactionTest.get_coefficient("GLC_b"))

    def test_check_mass_balance(self):
        # TODO: add cases
        pass

    def test_add_reactions_from_file(self):
        test_model = cb.Model(0)
        cr.add_reactions_from_file(
            model=test_model,
            filename=dir_input.joinpath("rxnToAdd_01_normal.txt"),
            directory=dir_data,
            database="META",
        )
        test_reactions = ["GLC_cb", "RXN_14462_p"]
        for test in test_reactions:
            self.assertTrue(
                test in [reaction.id for reaction in test_model.reactions]
            )

    def test_create_object(self):
        # CASE 1: metabolite from metacyc
        test_object = cr.create_object(
            identifier="GLC",
            directory=dir_data,
            compartment="c",
            database="META",
        )
        self.assertIsInstance(obj=test_object, cls=cb.Metabolite)
        self.assertEqual(first=test_object.compartment, second="c")
        # CASE 2: metabolite from kegg
        test_object = cr.create_object(
            identifier="C00026",
            directory=dir_data,
            compartment="p",
            database="KEGG",
        )
        self.assertIsInstance(obj=test_object, cls=cb.Metabolite)
        self.assertEqual(first=test_object.compartment, second="p")
        # CASE 3: reaction from metacyc
        test_object = cr.create_object(
            identifier="2.6.1.58-RXN",
            directory=dir_data,
            compartment="p",
            database="META",
        )
        self.assertIsInstance(obj=test_object, cls=cb.Reaction)
        self.assertIn(member="p", container=test_object.compartments)
        # CASE 4: reaction from KEGG
        test_object = cr.create_object(
            identifier="R02736",
            directory=dir_data,
            compartment="c",
            database="KEGG",
        )
        self.assertIsInstance(obj=test_object, cls=cb.Reaction)
        self.assertIn(member="c", container=test_object.compartments)


if __name__ == "__main__":
    unittest.main(verbosity=2)
