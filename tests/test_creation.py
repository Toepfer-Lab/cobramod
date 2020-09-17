#!/usr/bin/env python3
from pathlib import Path
from cobramod import creation as cr
import unittest
import xml.etree.ElementTree as ET
import cobra as cb

dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_biocyc = Path.cwd().joinpath("tests").joinpath("data").joinpath("biocyc")

if not dir_biocyc.exists():
    dir_biocyc.mkdir(parents=True)


class ModulTesting(unittest.TestCase):
    # TODO: use replacement dictionaries !!

    def test_get_xml_from_biocyc(self):
        # CASE 1: Directory does not exist
        self.assertRaises(
            FileNotFoundError,
            cr.get_xml_from_biocyc,
            # args
            directory=Path.cwd().joinpath("noDIr"),
            identifier="WATER",
            database="META"
        )
        # CASE 2: ID not found
        self.assertRaises(
            Warning,
            cr.get_xml_from_biocyc,
            directory=dir_biocyc,
            identifier="WATE",
            database="META"
        )
        # CASE 3: Proper usage with ET.Element
        dir_biocyc.joinpath("META").joinpath("WATER.xml").unlink()
        self.assertIsInstance(
            cr.get_xml_from_biocyc(
                directory=dir_biocyc,
                identifier="WATER",
                database="META"),
            ET.Element)
        # CASE 4: not found in database
        if dir_biocyc.joinpath("CPD-15326.xml").exists():
            dir_biocyc.joinpath("CPD-15326.xml").unlink()
        self.assertRaises(
            Warning,
            cr.get_xml_from_biocyc,
            directory=dir_biocyc,
            identifier="CPD-15326",
            database="ARA"
        )

    def test__create_meta_from_string(self):
        # CASE 1: Correct input
        testInput = "MALTOSE_b, MALTOSE[b], b, C12H22O11, 0"
        testMeta = cr._create_meta_from_string(line_string=testInput)
        # Checking that new meta is in Model
        self.assertIsInstance(testMeta, cb.Metabolite)
        # TODO: verify attributes

    def test_create_meta_from_root(self):
        # Case 0a: not string
        self.assertRaises(
            TypeError, cr.create_meta_from_root, root=list())
        # CASE 1: Correct input in cytosol
        testMeta = cr.create_meta_from_root(
            root="HOMOMETHIONINE", directory=dir_biocyc, database="META")
        self.assertIsInstance(testMeta, cb.Metabolite)
        # TODO: verify attributes

    def test_add_meta_from_string(self):
        # CASE 1: Not found in ARA, derivate to META
        test_metabolite = cr.add_meta_from_string(
            line="CPD-15326, c", model=cb.Model(0),
            database="ARA",
            directory=dir_biocyc)
        self.assertIsInstance(test_metabolite, cb.Metabolite)
        # CASE 2: Replace metabolite
        test_metabolite = cr.add_meta_from_string(
            line="Amino-Acids-20, p", model=cb.Model(0),
            database="ARA",
            directory=dir_biocyc, replacement_dict={
                "Amino-Acids-20": "GLT"})
        self.assertEqual(test_metabolite.id, "GLT_p")

    def test_add_meta_from_file(self):
        # FIXME: fix blank spaces
        # Creating copy
        testModel = cb.Model(0)
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
            testModel,
            Path.cwd().joinpath("nofile"),
        )
        # CASE 1: Metabolite is not found (or misspelled)
        self.assertRaises(
            Warning,
            cr.add_meta_from_file,
            model=testModel,
            filename=dir_input.joinpath("metaToAdd_02_misspelled.txt"),
            # Directory to save / check for xml files
            directory=dir_biocyc, database="META"
        )
        # CASE 2: Bad format (e. g. charge is missing).
        self.assertRaises(
            IndexError,
            cr.add_meta_from_file,
            model=testModel,
            filename=dir_input.joinpath("metaToAdd_03_badFormat.txt"),
            # Directory to save / check for xml files
            directory=dir_biocyc, database="META"
        )
        # CASE 3: Normal input
        cr.add_meta_from_file(
            model=testModel,
            # File with Metabolites to add
            filename=dir_input.joinpath("metaToAdd_01_normal.txt"),
            # Directory to save / check for xml files
            directory=dir_biocyc, database="META"
        )
        # Length should increase
        self.assertEqual(
            len(testModel.metabolites), 2)
        # Both Names in metabolites
        test_list = ["HOMOMETHIONINE_c", "MALTOSE_b"]
        self.assertTrue(
            all([x in [meta.id for meta in testModel.metabolites]
                for x in test_list]))

    def test__create_sides_for_reaction(self):
        # NOTE: Check for compartments !!
        # CASE 0a: from argument for 'side'
        self.assertRaises(
            Warning,
            cr._create_sides_for_reaction,
            model=cb.Model,
            root=ET.Element,
            temp_reaction=cb.Reaction,
            side="RiFj")
        # CASE 1a: left side a+b -> c+d
        test_root = cr.get_xml_from_biocyc(  # root
            directory=dir_biocyc, identifier="OXALODECARB-RXN",
            database="META")
        test_list = ["OXALACETIC_ACID_c", "PROTON_c"]
        test_reaction = cr._create_sides_for_reaction(
            model=cb.Model(0), root=test_root, temp_reaction=cb.Reaction(0),
            side="left", directory=dir_biocyc, database="META")
        # size should be 2
        self.assertEqual(2, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in test_list])
        )
        # CASE 1b: right side a+b -> c+d
        test_reaction = cr._create_sides_for_reaction(
            model=cb.Model(0), root=test_root, temp_reaction=cb.Reaction(0),
            side="right", directory=dir_biocyc, database="META")
        test_list = ["PYRUVATE_c", "CARBON_DIOXIDE_c"]
        self.assertEqual(2, len(test_reaction.metabolites))
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in test_list])
        )
        # CASE 2a: left side a+3b -> c+d+e+2f
        test_root = cr.get_xml_from_biocyc(  # root
            directory=dir_biocyc,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="META")
        test_list = ["WATER_c", "GTP_c"]
        test_reaction = cr._create_sides_for_reaction(
            model=cb.Model(0), root=test_root, temp_reaction=cb.Reaction(0),
            side="left", directory=dir_biocyc, database="META")
        # size should be 2
        self.assertEqual(2, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in test_list]))
        # WATER should have a coefficient of -3
        self.assertEqual(
            -3, test_reaction.get_coefficient("WATER_c")
        )
        # CASE 2b: right side of a+3b -> c+d+e+2f
        test_reaction = cr._create_sides_for_reaction(
            model=cb.Model(0), root=test_root, temp_reaction=cb.Reaction(0),
            side="right", directory=dir_biocyc, database="META")
        test_list = [
            "PROTON_c", "PPI_c",
            "DIAMINO_OH_PHOSPHORIBOSYLAMINO_PYR_c", "FORMATE_c"]
        # size should be 4
        self.assertEqual(4, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in test_list]))
        # CASE 3: root a+b(Protein)+c -> c+d(Protein)+c
        test_root = cr.get_xml_from_biocyc(  # root
            directory=dir_biocyc,
            identifier="RXN-11413",
            database="META")
        test_list = [
            "INDOLE_3_ACETALDOXIME_c", "Red_NADPH_Hemoprotein_Reductases_c",
            "OXYGEN_MOLECULE_c"]
        test_reaction = cr._create_sides_for_reaction(
            model=cb.Model(0), root=test_root, temp_reaction=cb.Reaction(0),
            side="left", directory=dir_biocyc, database="META")
        # size should be 2
        self.assertEqual(3, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in test_list]))

    def test_build_reaction_from_xml(self):
        # CASE 1:  Reaction not found in Ara, derive to META
        if dir_biocyc.joinpath(
                "GLUTAMINE--PYRUVATE-AMINOTRANSFERASE-RXN.xml").exists():
            dir_biocyc.joinpath(
                "GLUTAMINE--PYRUVATE-AMINOTRANSFERASE-RXN.xml").unlink()
        test_reaction = cr.build_reaction_from_xml(
            root="GLUTAMINE--PYRUVATE-AMINOTRANSFERASE-RXN",
            directory=dir_biocyc, model=cb.Model(0), database="META",
            compartment="c"
        )
        self.assertIsInstance(test_reaction, cb.Reaction)
        # CASE 2: Replacing name
        test_reaction = cr.build_reaction_from_xml(
            root="RXN-2205", directory=dir_biocyc, database="META",
            compartment="c",
            model=cb.Model(0), replacement_dict={
                "Amino-Acids-20": "GLT"})
        self.assertIn(
            "GLT_c", [meta.id for meta in test_reaction.metabolites])

    def test__add_reaction_line_to_model(self):
        # CASE 1: using delimiter, compartment is cytosol
        test_model = cb.Model(0)
        test_line = (
            'RXN_17742_c, RXN_17742_c |'
            'Oxidized-ferredoxins_c:-1, Reduced-ferredoxins_c: 1')
        cr._add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dir_biocyc,
            database="META")
        self.assertTrue(
            "RXN_17742_c" in [
                reaction.id for reaction in test_model.reactions])
        # CASE 2: No delimiter
        test_model = cb.Model(0)
        test_line = 'RXN-14462, c'
        cr._add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dir_biocyc,
            database="META")
        self.assertTrue(
            "RXN_14462_c" in [
                reaction.id for reaction in test_model.reactions])
        # CASE 2: No delimiter, compartment p
        test_model = cb.Model(0)
        test_line = 'RXN-14462, p'
        cr._add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dir_biocyc,
            database="META")
        self.assertTrue(
            "RXN_14462_p" in [
                reaction.id for reaction in test_model.reactions])

    def test_add_reaction_from_root(self):
        # CASE 0: Model invalid
        self.assertRaises(
            TypeError, cr.add_reaction_from_root, "NotModel", int(),
            database="META", directory=dir_biocyc)
        # CASE 1: root invalid
        self.assertRaises(
            TypeError, cr.add_reaction_from_root, cb.Model(0), int(),
            database="META", directory=dir_biocyc)
        # CASE 2a: proper root (ET.Element) a+b -> c+d
        # creating test root
        test_root = cr.get_xml_from_biocyc(
            directory=dir_biocyc,
            identifier="OXALODECARB-RXN", database="META")
        testModel = cb.Model(0)
        cr.add_reaction_from_root(
            model=testModel, root=test_root,
            directory=dir_biocyc, database="META")
        # Testing for new metabolites in Reaction
        self.assertTrue(
            "OXALODECARB_RXN_c" in [rxn.id for rxn in testModel.reactions])
        test_list = [
            "PYRUVATE_c", "CARBON_DIOXIDE_c", "OXALACETIC_ACID_c", "PROTON_c"]
        meta_list = [
            meta.id for meta in testModel.reactions.get_by_id(
                "OXALODECARB_RXN_c").metabolites]
        self.assertTrue(
            all([test in meta_list for test in test_list]))
        # CASE 2b: as with 2a but a string as root
        # Reaction
        testModel = cb.Model(0)
        cr.add_reaction_from_root(
            model=testModel, root="OXALODECARB-RXN",
            directory=dir_biocyc, database="META")
        # Testing for new metabolites in Reaction
        self.assertTrue(
            "OXALODECARB_RXN_c" in [rxn.id for rxn in testModel.reactions])
        test_list = [
            "PYRUVATE_c", "CARBON_DIOXIDE_c", "OXALACETIC_ACID_c", "PROTON_c"]
        meta_list = [
            meta.id for meta in testModel.reactions.get_by_id(
                "OXALODECARB_RXN_c").metabolites]
        self.assertTrue(
            all([test in meta_list for test in test_list]))

        # CASE 3: a+3b -> c+d+e+2f
        test_list = [
            "WATER_c", "GTP_c", "PROTON_c", "PPI_c",
            "DIAMINO_OH_PHOSPHORIBOSYLAMINO_PYR_c", "FORMATE_c"]
        testModel = cb.Model(0)
        cr.add_reaction_from_root(
            model=testModel, root="GTP-CYCLOHYDRO-II-RXN",
            directory=dir_biocyc, database="META")
        meta_list = [
            meta.id for meta in testModel.reactions.get_by_id(
                "GTP_CYCLOHYDRO_II_RXN_c").metabolites]
        self.assertTrue(
            "GTP_CYCLOHYDRO_II_RXN_c" in [
                rxn.id for rxn in testModel.reactions])
        # size should be 6
        self.assertEqual(
            6, len(testModel.reactions.get_by_id(
                "GTP_CYCLOHYDRO_II_RXN_c").metabolites))
        # names should be there
        self.assertTrue(
            all([test in meta_list for test in test_list]))
        # PROTON should have a coefficient of 2
        self.assertEqual(
            2, testModel.reactions.get_by_id(
                "GTP_CYCLOHYDRO_II_RXN_c").get_coefficient("PROTON_c"))
        # WATER should have a coefficient of -3
        self.assertEqual(
            -3, testModel.reactions.get_by_id(
                "GTP_CYCLOHYDRO_II_RXN_c").get_coefficient("WATER_c")
        )
        # FIXME: check reversibility
        # CASE 5: root a+b(Protein)+c -> d+e(Protein)+f
        testModel = cb.Model(0)
        test_list = [
            "CPD_12386_c", "Ox_NADPH_Hemoprotein_Reductases_c",
            "WATER_c", "INDOLE_3_ACETALDOXIME_c",
            "Red_NADPH_Hemoprotein_Reductases_c", "OXYGEN_MOLECULE_c"]
        cr.add_reaction_from_root(
            model=testModel, root="RXN-11413",
            directory=dir_biocyc, database="META")
        meta_list = [
            meta.id for meta in testModel.reactions.get_by_id(
                "RXN_11413_c").metabolites]
        self.assertTrue(
            "RXN_11413_c" in [
                rxn.id for rxn in testModel.reactions])
        # size should be 6
        self.assertEqual(6, len(testModel.reactions.get_by_id(
                "RXN_11413_c").metabolites))
        # names should be there
        self.assertTrue(
            all([test in meta_list for test in test_list]))
        # CASE 6: root a+b(Protein)+c -> d+e(Protein)+f (compartment p)
        testModel = cb.Model(0)
        test_list = [
            "CPD_12386_p", "Ox_NADPH_Hemoprotein_Reductases_p",
            "WATER_p", "INDOLE_3_ACETALDOXIME_p",
            "Red_NADPH_Hemoprotein_Reductases_p", "OXYGEN_MOLECULE_p"]
        cr.add_reaction_from_root(
            model=testModel, root="RXN-11413",
            directory=dir_biocyc, compartment="p", database="META")
        meta_list = [
            meta.id for meta in testModel.reactions.get_by_id(
                "RXN_11413_p").metabolites]
        self.assertTrue(
            "RXN_11413_p" in [
                rxn.id for rxn in testModel.reactions])
        # size should be 6
        self.assertEqual(6, len(testModel.reactions.get_by_id(
                "RXN_11413_p").metabolites))
        # names should be there
        self.assertTrue(
            all([test in meta_list for test in test_list]))
        # CASE 7: str not found
        self.assertRaises(
            Warning, cr.add_reaction_from_root, model=cb.Model(0),
            root="fakeroot", directory=dir_biocyc, database="META")

    def test__read_lines(self):
        # CASE 0: Comments and blank lines
        with open(
                file=dir_input.joinpath("test_reading_lines.txt")) as f:
            line = list(cr._read_lines(f=f))
        self.assertEqual(len(line), 4)

    def test__build_dict_for_metabolites(self):
        # CASE 0a: TypeError
        self.assertRaises(TypeError, cr._build_dict_for_metabolites, str())
        # CASE 0b: Coefficient missing
        self.assertRaises(
            ValueError, cr._build_dict_for_metabolites,
            string_list=[" GLC_c:"])
        self.assertDictEqual(
            {'GLC_b': 1.0, 'GLC_c': -1.0},
            cr._build_dict_for_metabolites(
                string_list=[" GLC_c:-1", "GLC_b: 1"]))

    def test_create_custom_reaction(self):
        # CASE 0: wrong format, no delimiter
        self.assertRaises(
            IndexError, cr.create_custom_reaction,
            line_string="GLC_cb, GLC_cb GLC_c:-1, GLC_b:1", model=cb.Model(0))
        # CASE 1: No ID detected
        self.assertRaises(
            Warning, cr.create_custom_reaction,
            line_string=" |GLC_c:-1, GLC_b:1", model=cb.Model(0))
        # CASE 2: Normal, ID and name differ
        ReactionTest = cr.create_custom_reaction(
            line_string="GLC_cb, Glucose Transport|GLC_c:-1, GLC_b:1",
            model=cb.Model(0), directory=dir_biocyc, database="META")
        # Checking if ID, name and instance are correct.
        self.assertTrue(
            all([
                ReactionTest.id == "GLC_cb",
                ReactionTest.name == "Glucose Transport",
                isinstance(ReactionTest, cb.Reaction)])
        )
        test_list = ["GLC_c", "GLC_b"]
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in ReactionTest.metabolites]
                for x in test_list]))
        self.assertEqual(
            -1, ReactionTest.get_coefficient("GLC_c"))
        self.assertEqual(
            1, ReactionTest.get_coefficient("GLC_b"))

    def test_check_mass_balance(self):
        # TODO: !!
        pass

    def test_add_reactions_from_file(self):
        testModel = cb.Model(0)
        cr.add_reactions_from_file(
            model=testModel,
            filename=dir_input.joinpath("rxnToAdd_01_normal.txt"),
            directory=dir_biocyc, database="META")
        test_reactions = ["GLC_cb", "RXN_14462_p"]
        for test in test_reactions:
            self.assertTrue(
                test in [reaction.id for reaction in testModel.reactions])


if __name__ == "__main__":
    unittest.main(verbosity=2)
