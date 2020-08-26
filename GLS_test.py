#!/usr/bin/env python3
from pathlib import Path
import GLS
import unittest
import xml.etree.ElementTree as ET
import cobra as cb

dirBiocyc = Path.cwd().joinpath("biocyc")
dirInput = Path.cwd().joinpath("input")


class ModulTesting(unittest.TestCase):
    # TODO: create tests for more complex models!!
    # TODO: use replacement dictionaries !!
    # TODO: check for subdatabase in Biocyc

    def test_get_xml_from_biocyc(self):
        # CASE 1: Directory does not exist
        self.assertRaises(
            FileNotFoundError,
            GLS.get_xml_from_biocyc,
            # args
            directory=Path.cwd().joinpath("noDIr"),
            bioID="WATER",
        )
        # CASE 2: ID not found
        self.assertRaises(
            Warning,
            GLS.get_xml_from_biocyc,
            directory=dirBiocyc,
            bioID="WATE",
        )
        # CASE 3: Proper usage with ET.Element
        dirBiocyc.joinpath("WATER.xml").unlink()
        self.assertIsInstance(
            GLS.get_xml_from_biocyc(
                directory=dirBiocyc,
                bioID="WATER"),
            ET.Element)

    def test_create_meta_from_string(self):
        # CASE 1: Correct input
        testInput = "MALTOSE_b, MALTOSE[b], b, C12H22O11, 0"
        testMeta = GLS.create_meta_from_string(line_string=testInput)
        # Checking that new meta is in Model
        self.assertIsInstance(testMeta, cb.Metabolite)
        # TODO: verify attributes

    def test_create_meta_from_root(self):
        # Case 0a: not string
        self.assertRaises(
            TypeError, GLS.create_meta_from_root, root=list())
        # CASE 1: Correct input in cytosol
        testMeta = GLS.create_meta_from_root(
            "HOMOMETHIONINE", directory=dirBiocyc)
        self.assertIsInstance(testMeta, cb.Metabolite)
        # TODO: verify attributes

    def test_add_meta_from_file(self):
        # FIXME: fix blank spaces
        # Creating copy
        testModel = cb.Model(0)
        # Testing if model is not cobra.Model
        self.assertRaises(
            TypeError,
            GLS.add_meta_from_file,
            "notModel",
            Path.cwd().joinpath("nofile"),
        )
        # Testing if file is not found
        self.assertRaises(
            FileNotFoundError,
            GLS.add_meta_from_file,
            testModel,
            Path.cwd().joinpath("nofile"),
        )
        # CASE 1: Metabolite is not found (or misspelled)
        self.assertRaises(
            Warning,
            GLS.add_meta_from_file,
            model=testModel,
            filename=dirInput.joinpath("metaToAdd_02_misspelled.txt"),
            # Directory to save / check for xml files
            directory=dirBiocyc
        )
        # CASE 2: Bad format (e. g. charge is missing).
        self.assertRaises(
            IndexError,
            GLS.add_meta_from_file,
            model=testModel,
            filename=dirInput.joinpath("metaToAdd_03_badFormat.txt"),
            # Directory to save / check for xml files
            directory=dirBiocyc
        )
        # CASE 3: Normal input
        GLS.add_meta_from_file(
            model=testModel,
            # File with Metabolites to add
            filename=dirInput.joinpath("metaToAdd_01_normal.txt"),
            # Directory to save / check for xml files
            directory=dirBiocyc
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
            GLS.create_sides_for_reaction,
            model=cb.Model,
            root=ET.Element,
            temp_reaction=cb.Reaction,
            side="RiFj")
        # CASE 1a: left side a+b -> c+d 
        testRoot = GLS.get_xml_from_biocyc(  # root
            dirBiocyc,
            "OXALODECARB-RXN")
        toCheckList = ["OXALACETIC_ACID_c", "PROTON_c"]
        test_reaction = GLS.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="left", directory=dirBiocyc) # !!
        # size should be 2
        self.assertEqual(2, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList])
        )
        # CASE 1b: right side a+b -> c+d 
        test_reaction = GLS.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="right", directory=dirBiocyc)
        toCheckList = ["PYRUVATE_c", "CARBON_DIOXIDE_c"]
        self.assertEqual(2, len(test_reaction.metabolites))
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList])
        )
        # CASE 2a: left side a+3b -> c+d+e+2f
        testRoot = GLS.get_xml_from_biocyc(  # root
            dirBiocyc,
            "GTP-CYCLOHYDRO-II-RXN")
        toCheckList = ["WATER_c", "GTP_c"]
        test_reaction = GLS.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="left", directory=dirBiocyc)
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
        test_reaction = GLS.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="right", directory=dirBiocyc)
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
        testRoot = GLS.get_xml_from_biocyc(  # root
            dirBiocyc,
            "RXN-11413")
        toCheckList = [
            "INDOLE_3_ACETALDOXIME_c", "Red_NADPH_Hemoprotein_Reductases_c",
            "OXYGEN_MOLECULE_c"]
        test_reaction = GLS.create_sides_for_reaction(
            model=cb.Model(0), root=testRoot, temp_reaction=cb.Reaction(0),
            side="left", directory=dirBiocyc)
        # size should be 2
        self.assertEqual(3, len(test_reaction.metabolites))
        # names should be there
        self.assertTrue(
            all([x in [meta.id for meta in test_reaction.metabolites]
                for x in toCheckList]))
    def test_build_reaction_from_xml(self):
        # NOTE: finish this part
        test_model = cb.Model(0)
        

    def test_add_reaction_line_to_model(self):
        # CASE 1: using delimiter, compartment is cytosol
        test_model = cb.Model(0)
        test_line = (
            f'RXN_17742_c, RXN_17742_c |'
            f'Oxidized-ferredoxins_c:-1, Reduced-ferredoxins_c: 1')
        GLS.add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dirBiocyc)
        self.assertTrue(
            "RXN_17742_c" in [
                reaction.id for reaction in test_model.reactions])
        # CASE 2: No delimiter
        test_model = cb.Model(0)
        test_line = f'RXN-14462, c'
        GLS.add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dirBiocyc)
        self.assertTrue(
            "RXN_14462_c" in [
                reaction.id for reaction in test_model.reactions])
        # CASE 2: No delimiter, compartment p
        test_model = cb.Model(0)
        test_line = f'RXN-14462, p'
        GLS.add_reaction_line_to_model(
            line=test_line, model=test_model, directory=dirBiocyc)
        self.assertTrue(
            "RXN_14462_p" in [
                reaction.id for reaction in test_model.reactions])
    
    def test_add_reaction_from_root(self):
        # CASE 0: Model invalid
        self.assertRaises(
            TypeError, GLS.add_reaction_from_root, "NotModel", int())
        # CASE 1: root invalid
        self.assertRaises(
            TypeError, GLS.add_reaction_from_root, cb.Model(0), int())
        # CASE 2a: proper root (ET.Element) a+b -> c+d
        # creating test root
        testRoot = GLS.get_xml_from_biocyc(
            directory=dirBiocyc,
            bioID="OXALODECARB-RXN")
        testModel = cb.Model(0)
        GLS.add_reaction_from_root(
            model=testModel, root=testRoot,
            directory=dirBiocyc)
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
        GLS.add_reaction_from_root(
            model=testModel, root="OXALODECARB-RXN",
            directory=dirBiocyc)
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
        GLS.add_reaction_from_root(
            model=testModel, root="GTP-CYCLOHYDRO-II-RXN",
            directory=dirBiocyc)
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
        GLS.add_reaction_from_root(
            model=testModel, root="RXN-11413",
            directory=dirBiocyc)
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
        GLS.add_reaction_from_root(
            model=testModel, root="RXN-11413",
            directory=dirBiocyc, compartment="p")
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
            Warning, GLS.add_reaction_from_root, model=cb.Model(0),
            root="fakeroot", directory=dirBiocyc)

    def test_read_lines(self):
        # CASE 0: Comments and blank lines
        with open(
            file=dirInput.joinpath("test_reading_lines.txt")) as f:
            line = list(GLS.read_lines(f=f))
        self.assertEqual(len(line), 4)

    def test_createDictForMetabolites(self):
        # CASE 0a: TypeError
        self.assertRaises(TypeError, GLS.createDictForMetabolites, str())
        # CASE 0b: Coefficient missing
        self.assertRaises(
            ValueError, GLS.createDictForMetabolites,
            metaString=[" GLC_c:"])
        self.assertDictEqual(
            {'GLC_b': 1.0, 'GLC_c': -1.0},
            GLS.createDictForMetabolites(
                metaString=[" GLC_c:-1", "GLC_b: 1"]))

    def test_create_custom_reaction(self):
        # CASE 0: wrong format, no delimiter
        self.assertRaises(
            IndexError, GLS.create_custom_reaction,
            line_string="GLC_cb, GLC_cb GLC_c:-1, GLC_b:1")
        # CASE 1: No ID detected
        self.assertRaises(
            Warning, GLS.create_custom_reaction,
            line_string="    |GLC_c:-1, GLC_b:1")
        # CASE 2: Normal, ID and name differ
        ReactionTest = GLS.create_custom_reaction(
            line_string="GLC_cb, Glucose Transport|GLC_c:-1, GLC_b:1",
            model=cb.Model(0), directory=dirBiocyc)
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
        GLS.add_reaction_from_file(
            model=testModel,
            filename=dirInput.joinpath("rxnToAdd_01_normal.txt"),
            directory=dirBiocyc)
        test_reactions = ["GLC_cb", "RXN_14462_p"]
        for test in test_reactions:
            self.assertTrue(
                test in [reaction.id for reaction in testModel.reactions])


if __name__ == "__main__":
    unittest.main(verbosity=2)
