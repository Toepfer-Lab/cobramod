#!/usr/bin/env python3
from pathlib import Path
import creation as cr
import unittest
import xml.etree.ElementTree as ET
import cobra as cb

dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_biocyc = Path.cwd().joinpath("tests").joinpath("data").joinpath("biocyc")

if not dir_biocyc.exists():
    dir_biocyc.mkdir(parents=True)


class ModulTesting(unittest.TestCase):
    # TODO: create tests for more complex models!!
    # TODO: use replacement dictionaries !!
    # TODO: check for subdatabase in Biocyc

    def test_get_xml_from_biocyc(self):
        # CASE 1: Directory does not exist
        self.assertRaises(
            FileNotFoundError,
            cr.get_xml_from_biocyc,
            # args
            directory=Path.cwd().joinpath("noDIr"),
            bioID="WATER",
        )
        # CASE 2: ID not found
        self.assertRaises(
            Warning,
            cr.get_xml_from_biocyc,
            directory=dir_biocyc,
            bioID="WATE",
        )
        # CASE 3: Proper usage with ET.Element
        dir_biocyc.joinpath("WATER.xml").unlink()
        self.assertIsInstance(
            cr.get_xml_from_biocyc(
                directory=dir_biocyc,
                bioID="WATER"),
            ET.Element)

    def test_create_meta_from_string(self):
        # CASE 1: Correct input
        testInput = "MALTOSE_b, MALTOSE[b], b, C12H22O11, 0"
        testMeta = cr.create_meta_from_string(line_string=testInput)
        # Checking that new meta is in Model
        self.assertIsInstance(testMeta, cb.Metabolite)
        # TODO: verify attributes

    def test_create_meta_from_root(self):
        # Case 0a: not string
        self.assertRaises(
            TypeError, cr.create_meta_from_root, root=list())
        # CASE 1: Correct input in cytosol
        testMeta = cr.create_meta_from_root(
            "HOMOMETHIONINE", directory=dir_biocyc)
        self.assertIsInstance(testMeta, cb.Metabolite)
        # TODO: verify attributes

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
            directory=dir_biocyc
        )
        # CASE 2: Bad format (e. g. charge is missing).
        self.assertRaises(
            IndexError,
            cr.add_meta_from_file,
            model=testModel,
            filename=dir_input.joinpath("metaToAdd_03_badFormat.txt"),
            # Directory to save / check for xml files
            directory=dir_biocyc
        )
        # CASE 3: Normal input
        cr.add_meta_from_file(
            model=testModel,
            # File with Metabolites to add
            filename=dir_input.joinpath("metaToAdd_01_normal.txt"),
            # Directory to save / check for xml files
            directory=dir_biocyc
        )
        # Length should increase
        self.assertEqual(
            len(testModel.metabolites), 2)
        # Both Names in metabolites
        toCheckList = ["HOMOMETHIONINE_c", "MALTOSE_b"]
        self.assertTrue(
            all([x in [meta.id for meta in testModel.metabolites]
                for x in toCheckList]))

    def test_create_sides_for_reaction(self):
        # NOTE: Check for compartments !!
        # CASE 0a: from argument for 'side'
        self.assertRaises(
            Warning,
            cr.create_sides_for_reaction,
            model=cb.Model,
            root=ET.Element,
            temp_reaction=cb.Reaction,
            side="RiFj")
        # CASE 1a: left side a+b -> c+d
        testRoot = cr.get_xml_from_biocyc(  # root
            dir_biocyc,
            "OXALODECARB-RXN")
        toCheckList = ["OXALACETIC_ACID_c", "PROTON_c"]
        test_reaction = cr.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="left", directory=dir_biocyc)
        # size should be 2
        self.assertEqual(2, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList])
        )
        # CASE 1b: right side a+b -> c+d
        test_reaction = cr.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="right", directory=dir_biocyc)
        toCheckList = ["PYRUVATE_c", "CARBON_DIOXIDE_c"]
        self.assertEqual(2, len(test_reaction.metabolites))
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList])
        )
        # CASE 2a: left side a+3b -> c+d+e+2f
        testRoot = cr.get_xml_from_biocyc(  # root
            dir_biocyc,
            "GTP-CYCLOHYDRO-II-RXN")
        toCheckList = ["WATER_c", "GTP_c"]
        test_reaction = cr.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="left", directory=dir_biocyc)
        # size should be 2
        self.assertEqual(2, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList]))
        # WATER should have a coefficient of -3
        self.assertEqual(
            -3, test_reaction.get_coefficient("WATER_c")
        )
        # CASE 2b: right side of a+3b -> c+d+e+2f
        test_reaction = cr.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="right", directory=dir_biocyc)
        toCheckList = [
            "PROTON_c", "PPI_c",
            "DIAMINO_OH_PHOSPHORIBOSYLAMINO_PYR_c", "FORMATE_c"]
        # size should be 4
        self.assertEqual(4, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList]))
        # CASE 3: root a+b(Protein)+c -> c+d(Protein)+c
        testRoot = cr.get_xml_from_biocyc(  # root
            dir_biocyc,
            "RXN-11413")
        toCheckList = [
            "INDOLE_3_ACETALDOXIME_c", "Red_NADPH_Hemoprotein_Reductases_c",
            "OXYGEN_MOLECULE_c"]
        test_reaction = cr.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="left", directory=dir_biocyc)
        # size should be 2
        self.assertEqual(3, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList]))

    def test_build_reaction_from_xml(self):
        # TODO: finish this part
        # test_model = cb.Model(0)
        pass

    def test_add_reaction_line_to_model(self):
        # CASE 1: using delimiter, compartment is cytosol
        test_model = cb.Model(0)
        test_line = (
            'RXN_17742_c, RXN_17742_c |'
            'Oxidized-ferredoxins_c:-1, Reduced-ferredoxins_c: 1')
        cr.add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dir_biocyc)
        self.assertTrue(
            "RXN_17742_c" in [
                reaction.id for reaction in test_model.reactions])
        # CASE 2: No delimiter
        test_model = cb.Model(0)
        test_line = 'RXN-14462, c'
        cr.add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dir_biocyc)
        self.assertTrue(
            "RXN_14462_c" in [
                reaction.id for reaction in test_model.reactions])
        # CASE 2: No delimiter, compartment p
        test_model = cb.Model(0)
        test_line = 'RXN-14462, p'
        cr.add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dir_biocyc)
        self.assertTrue(
            "RXN_14462_p" in [
                reaction.id for reaction in test_model.reactions])

    def test_add_reaction_from_root(self):
        # CASE 0: Model invalid
        self.assertRaises(
            TypeError, cr.add_reaction_from_root, "NotModel", int())
        # CASE 1: root invalid
        self.assertRaises(
            TypeError, cr.add_reaction_from_root, cb.Model(0), int())
        # CASE 2a: proper root (ET.Element) a+b -> c+d
        # creating test root
        testRoot = cr.get_xml_from_biocyc(
            directory=dir_biocyc,
            bioID="OXALODECARB-RXN")
        testModel = cb.Model(0)
        cr.add_reaction_from_root(
            model=testModel, root=testRoot,
            directory=dir_biocyc)
        # Testing for new metabolites in Reaction
        self.assertTrue(
            "OXALODECARB_RXN_c" in [rxn.id for rxn in testModel.reactions])
        toCheckList = [
            "PYRUVATE_c", "CARBON_DIOXIDE_c", "OXALACETIC_ACID_c", "PROTON_c"]
        meta_list = [
            meta.id for meta in testModel.reactions.get_by_id(
                "OXALODECARB_RXN_c").metabolites]
        self.assertTrue(
            all([test in meta_list for test in toCheckList]))
        # CASE 2b: as with 2a but a string as root
        # Reaction
        testModel = cb.Model(0)
        cr.add_reaction_from_root(
            model=testModel, root="OXALODECARB-RXN",
            directory=dir_biocyc)
        # Testing for new metabolites in Reaction
        self.assertTrue(
            "OXALODECARB_RXN_c" in [rxn.id for rxn in testModel.reactions])
        toCheckList = [
            "PYRUVATE_c", "CARBON_DIOXIDE_c", "OXALACETIC_ACID_c", "PROTON_c"]
        meta_list = [
            meta.id for meta in testModel.reactions.get_by_id(
                "OXALODECARB_RXN_c").metabolites]
        self.assertTrue(
            all([test in meta_list for test in toCheckList]))

        # CASE 3: a+3b -> c+d+e+2f
        toCheckList = [
            "WATER_c", "GTP_c", "PROTON_c", "PPI_c",
            "DIAMINO_OH_PHOSPHORIBOSYLAMINO_PYR_c", "FORMATE_c"]
        testModel = cb.Model(0)
        cr.add_reaction_from_root(
            model=testModel, root="GTP-CYCLOHYDRO-II-RXN",
            directory=dir_biocyc)
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
            all([test in meta_list for test in toCheckList]))
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
        toCheckList = [
            "CPD_12386_c", "Ox_NADPH_Hemoprotein_Reductases_c",
            "WATER_c", "INDOLE_3_ACETALDOXIME_c",
            "Red_NADPH_Hemoprotein_Reductases_c", "OXYGEN_MOLECULE_c"]
        cr.add_reaction_from_root(
            model=testModel, root="RXN-11413",
            directory=dir_biocyc)
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
            all([test in meta_list for test in toCheckList]))
        # CASE 6: root a+b(Protein)+c -> d+e(Protein)+f (compartment p)
        testModel = cb.Model(0)
        toCheckList = [
            "CPD_12386_p", "Ox_NADPH_Hemoprotein_Reductases_p",
            "WATER_p", "INDOLE_3_ACETALDOXIME_p",
            "Red_NADPH_Hemoprotein_Reductases_p", "OXYGEN_MOLECULE_p"]
        cr.add_reaction_from_root(
            model=testModel, root="RXN-11413",
            directory=dir_biocyc, compartment="p")
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
            all([test in meta_list for test in toCheckList]))
        # CASE 7: str not found
        self.assertRaises(
            Warning, cr.add_reaction_from_root, model=cb.Model(0),
            root="fakeroot", directory=dir_biocyc)

    def test_read_lines(self):
        # CASE 0: Comments and blank lines
        with open(
                file=dir_input.joinpath("test_reading_lines.txt")) as f:
            line = list(cr.read_lines(f=f))
        self.assertEqual(len(line), 4)

    def test_createDictForMetabolites(self):
        # CASE 0a: TypeError
        self.assertRaises(TypeError, cr.createDictForMetabolites, str())
        # CASE 0b: Coefficient missing
        self.assertRaises(
            ValueError, cr.createDictForMetabolites,
            metaString=[" GLC_c:"])
        self.assertDictEqual(
            {'GLC_b': 1.0, 'GLC_c': -1.0},
            cr.createDictForMetabolites(
                metaString=[" GLC_c:-1", "GLC_b: 1"]))

    def test_create_custom_reaction(self):
        # CASE 0: wrong format, no delimiter
        self.assertRaises(
            IndexError, cr.create_custom_reaction,
            line_string="GLC_cb, GLC_cb GLC_c:-1, GLC_b:1")
        # CASE 1: No ID detected
        self.assertRaises(
            Warning, cr.create_custom_reaction,
            line_string="    |GLC_c:-1, GLC_b:1")
        # CASE 2: Normal, ID and name differ
        ReactionTest = cr.create_custom_reaction(
            line_string="GLC_cb, Glucose Transport|GLC_c:-1, GLC_b:1",
            model=cb.Model(0), directory=dir_biocyc)
        # Checking if ID, name and instance are correct.
        self.assertTrue(
            all([
                ReactionTest.id == "GLC_cb",
                ReactionTest.name == "Glucose Transport",
                isinstance(ReactionTest, cb.Reaction)])
        )
        toCheckList = ["GLC_c", "GLC_b"]
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in ReactionTest.metabolites]
                for x in toCheckList]))
        self.assertEqual(
            -1, ReactionTest.get_coefficient("GLC_c"))
        self.assertEqual(
            1, ReactionTest.get_coefficient("GLC_b"))

    def test_stopAndShowMassBalance(self):
        # TODO: !!
        pass

    def test_add_reaction_from_file(self):
        testModel = cb.Model(0)
        cr.add_reaction_from_file(
            model=testModel,
            filename=dir_input.joinpath("rxnToAdd_01_normal.txt"),
            directory=dir_biocyc)
        test_reactions = ["GLC_cb", "RXN_14462_p"]
        for test in test_reactions:
            self.assertTrue(
                test in [reaction.id for reaction in testModel.reactions])


if __name__ == "__main__":
    unittest.main(verbosity=2)
