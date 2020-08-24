#!/usr/bin/env python3
from pathlib import Path
import unittest
import pathways
import xml.etree.ElementTree as ET
from GLS import get_xml_from_biocyc, add_reaction_from_root,\
    add_meta_line_to_model, add_meta_from_file, add_reaction_line_to_model
import cobra as cb

dirBiocyc = Path.cwd().joinpath("biocyc")
dirInput = Path.cwd().joinpath("input")


class ModulTesting(unittest.TestCase):
    @unittest.skip
    def test_dictRxnForPathway(self):
        self.assertRaises(
            TypeError, pathways.dictRxnForPathway, xmlRoot=str())
        self.assertIsInstance(
            pathways.dictRxnForPathway(
                get_xml_from_biocyc(bioID="PWY-1187", directory=dirBiocyc)),
            dict)
        # TODO: check for correct dictionary
    
    @unittest.skip
    def test_find_start_vertex(self):
        # CASE 1: regular lineal
        test_dict = {
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
            "D": ["E"]}
        test_result = list(pathways.find_start_vertex(
            vertex_dict=test_dict))
        # CASE 2: lineal with 2 paths
        test_dict = {
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
            "D": ["E"],
            "F": ["G"],
            "G": ["H", "Z"]
            }
        test_result = list(pathways.find_start_vertex(
            vertex_dict=test_dict))
        pass

    def test_find_end_vertex(self):
        # CASE 1: regular lineal
        test_dict = {
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
            "D": ["E"]}
        test_result = list(pathways.find_end_vertex(
            vertex_dict=test_dict))
        # CASE 2: lineal with 2 paths
        test_dict = {
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
            "D": ["E"],
            "F": ["G"],
            "G": ["H", "Z"]
            }
        test_result = list(pathways.find_end_vertex(
            vertex_dict=test_dict))
        pass

    @unittest.skip
    def test_get_edges_from_dict(self):
        # FIXME: find new solution for cyclical paths
        # CASE 1: regular
        start, ends = pathways.get_edges_from_dict({
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
            "D": ["E"]})
        self.assertTrue(all([
            len(start) == 1,
            len(ends) == 1]))
        # CASE 2: two ends
        start, ends = pathways.get_edges_from_dict({
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
            "D": ["E"],
            "G": ["E"]})
        self.assertTrue(all([
            len(start) == 1,  # E
            len(ends) == 2]))  # G, A
        # CASE 3: cycles
        start, ends = pathways.get_edges_from_dict({
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
            "D": ["E"],
            "E": ["A"]}
        )
        self.assertTrue(all([
            len(start) == 1,  # E
            len(ends) == 1,
            start == ends[0]]))
        # FIXME: add test
        # CASE 4: Cycle with 2 paths
        start, ends = pathways.get_edges_from_dict({
            "A": ["B"],
            "B": ["C"],
            "C": ["D", "G"],
            "D": ["E"],
            "E": ["A"],
            # In the middle
            "F": ["A"],
            "G": ["F"],
        })
        pass
    @unittest.skip
    def test_give_path_graph(self):
        # FIXME: add test
        # CASE 0: Lineal
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E"}
        test = pathways.give_path_graph(vertex="A", graph_dict=test_dict)
        # CASE 1: Lineal with two ends
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E", "G": "E"}
        start, ends = pathways.get_edges_from_dict(rootDict=test_dict)
        test_ends = [
            list(pathways.give_path_graph(
                vertex=vertex, graph_dict=test_dict)) for vertex in ends]
        # CASE 2: Cyclic
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "A",
            # In the middle
            "F": "A",
            "G": "F",
            "C": "G"
            }
        test_ends = list(pathways.give_path_graph(
                vertex="A", graph_dict=test_dict))
        # CASE 3: Cyclic with multiple ends

        pass
    @unittest.skip
    def test_createPathway(self):
        # CASE 0a: nothing passed
        self.assertRaises(TypeError, pathways.createPathway)
        # CASE 0b: both arguments passed
        self.assertRaises(
            Warning, pathways.createPathway,
            xmlRoot=ET.Element(0), customDict=dict())
        # CASE 1: Custom dictionary
        test_dict = {"A": "B", "B": "C", "C": "D", "D": "E", "G": "E"}
        listTest = pathways.createPathway(customDict=test_dict)
        toTest = [
            ["E", "D", "C", "B", "A"],
            ["E", "G"]]
        self.assertTrue(
            all([test in listTest for test in toTest])
        )
        # CASE 2: regular BioCycID
        listTest = pathways.createPathway(
            xmlRoot=get_xml_from_biocyc(bioID="PWY-1187", directory=dirBiocyc))
        # 3 end-metabolites + 1 in-between path (cyclic)
        self.assertEqual(len(listTest), 4)
        # CASE 3: cycle
        listTest = pathways.createPathway(
            xmlRoot=get_xml_from_biocyc(bioID="PWY-5690", directory=dirBiocyc))
        # Only one
        self.assertEqual(len(listTest[0]), 10)
        listTest = pathways.createPathway(
            xmlRoot=get_xml_from_biocyc(bioID="CALVIN-PWY", directory=dirBiocyc))
        pass

    @unittest.skip
    def test_create_reactions_for_list(self):
        # CASE 1: Regular pathway
        testList = pathways.createPathway(
            xmlRoot=get_xml_from_biocyc(
                bioID="PWY-1187", directory=dirBiocyc)
        )[0]  # only 1 path
        testRxns = pathways.create_reactions_for_list(
            pathwayList=testList, directory=dirBiocyc, model=cb.Model(0))
        self.assertTrue(all(
            [isinstance(rxn, cb.Reaction) for rxn in testRxns]
        ))
    @unittest.skip
    def test_find_next_demand(self):
        # CASE 0a: invalid Model
        self.assertRaises(
            TypeError, pathways.find_next_demand, model="Model",
            toCheckRxnID=str())
        # CASe 0b: no metabolites found
        testModel = cb.Model(0)
        reactionTest = cb.Reaction(id="test")
        testModel.add_reactions([reactionTest])
        self.assertRaises(
            Warning, pathways.find_next_demand, model=testModel,
            toCheckRxnID="test")
        # creating reaction and model
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="OXALODECARB-RXN",
            directory=dirBiocyc)
        # testModel.add_reactions([reactionTest])
        # CASE 1a: left --> right
        self.assertEqual(
            "PYRUVATE_c",
            pathways.find_next_demand(
                model=testModel, toCheckRxnID="OXALODECARB_RXN_c")
        )
        # CASE 1b: with a ignoreList
        self.assertEqual(
            "CARBON_DIOXIDE_c",
            pathways.find_next_demand(
                model=testModel, toCheckRxnID="OXALODECARB_RXN_c",
                ignoreList=["PYRUVATE_c"])
        )
        # CASE 2a: left <-- right
        testModel.reactions.get_by_id("OXALODECARB_RXN_c").bounds = (-1000, 0)
        self.assertEqual(
            "OXALACETIC_ACID_c",
            pathways.find_next_demand(
                model=testModel, toCheckRxnID="OXALODECARB_RXN_c")
        )
        # CASE 2b: with a ignoreList
        self.assertEqual(
            "PROTON_c",
            pathways.find_next_demand(
                model=testModel, toCheckRxnID="OXALODECARB_RXN_c",
                ignoreList=["OXALACETIC_ACID_c", "FAKEID"])
        )
        # NOTE: Test for reversible reactions
    @unittest.skip
    def test_isSink(self):
        # CASE 1: two normal reactions
        testModel = cb.Model(0)
        rxnList = pathways.create_reactions_for_list(
            model=testModel,
            pathwayList=["RXN-2206", "RXN-11414"],
            directory=dirBiocyc)
        testModel.add_reactions(rxnList)
        [testModel.add_boundary(
            meta, "sink") for meta in testModel.metabolites]
        toTest = ["RXN_2206_c", "RXN_11414_c"]
        self.assertTrue(not all([
            pathways.isSink(model=testModel, rxnID=reaction)
            for reaction in toTest]
        ))
        # CASE 2: two sinks
        toTest = [
            "SK_HOMOMETHIONINE_c", "SK_WATER_c"]
        self.assertTrue(all([
            pathways.isSink(model=testModel, rxnID=reaction)
            for reaction in toTest]
        ))
    @unittest.skip
    def test_fix_meta_for_boundaries(self):
        # CASE 0, ignore List
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="OXALODECARB-RXN",
            directory=dirBiocyc)
        self.assertRaises(
            Warning, pathways.fix_meta_for_boundaries, model=testModel,
            metabolite="PROTON_c", ignoreList=["PROTON_c"])
        # CASE 1: PROTON, two normal reactions, one after another
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="OXALODECARB-RXN",
            directory=dirBiocyc)
        pathways.fix_meta_for_boundaries(
            model=testModel, metabolite="PROTON_c")
        self.assertIn(
            "SK_PROTON_c", [
                reaction.id for reaction in testModel.metabolites.get_by_id(
                    "PROTON_c").reactions])
        add_reaction_from_root(
            model=testModel, xmlRoot="AROMATIC-L-AMINO-ACID-DECARBOXYLASE-RXN",
            directory=dirBiocyc)
        pathways.fix_meta_for_boundaries(
            model=testModel, metabolite="PROTON_c")
        self.assertNotIn(
            "SK_PROTON_c", [
                reaction.id for reaction in testModel.metabolites.get_by_id(
                    "PROTON_c").reactions])
        # CASE 2: Normal reaction, plus demand, plus test for sink
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="OXALODECARB-RXN",
            directory=dirBiocyc)
        testModel.add_boundary(
            testModel.metabolites.get_by_id("PROTON_c"), "demand")
        pathways.fix_meta_for_boundaries(
            model=testModel, metabolite="PROTON_c")
        self.assertIn(
            "DM_PROTON_c", [
                reaction.id for reaction in testModel.metabolites.get_by_id(
                    "PROTON_c").reactions])
        # CASE 3: Adding an extra reaction, nothing should be left but
        # reactions
        add_reaction_from_root(
            model=testModel, xmlRoot="AROMATIC-L-AMINO-ACID-DECARBOXYLASE-RXN",
            directory=dirBiocyc)
        pathways.fix_meta_for_boundaries(
            model=testModel, metabolite="PROTON_c")
        self.assertNotIn(
            "SK_PROTON_c", [
                reaction.id for reaction in testModel.metabolites.get_by_id(
                    "PROTON_c").reactions])
        self.assertNotIn(
            "DM_PROTON_c", [
                reaction.id for reaction in testModel.metabolites.get_by_id(
                    "PROTON_c").reactions])
    @unittest.skip
    def test_tryCreateOrRemoveSinkForSides(self):
        # CASE 0: invalid Model
        self.assertRaises(
            TypeError, pathways.tryCreateOrRemoveSinkForSides, model=str(),
            rxnID=str())
        self.assertRaises(
            ValueError, pathways.tryCreateOrRemoveSinkForSides,
            model=cb.Model(0),
            rxnID=str(), side="Not")
        # CASE 1: normal creation (left side)
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="OXALODECARB-RXN", directory=dirBiocyc)
        pathways.tryCreateOrRemoveSinkForSides(
            model=testModel, rxnID="OXALODECARB_RXN_c", side="left",
            ignoreList=["OXALACETIC_ACID_c"])
        self.assertTrue(
            "SK_PROTON_c" in [sink.id for sink in testModel.sinks]
        )
        # # CASE 2: Already sink
        pathways.tryCreateOrRemoveSinkForSides(
            model=testModel, rxnID="OXALODECARB_RXN_c", side="left",
            ignoreList=["OXALACETIC_ACID_c"])
        # # CASE 3: normal creation (right side)
        pathways.tryCreateOrRemoveSinkForSides(
            model=testModel, rxnID="OXALODECARB_RXN_c", side="right")
        self.assertTrue(
            "SK_PYRUVATE_c" in [sink.id for sink in testModel.sinks]
        )
    @unittest.skip
    def test_createAndCheckSinksForRxn(self):
        # CASE 0: Model invalid
        self.assertRaises(
            TypeError, pathways.createAndCheckSinksForRxn, model=str(),
            rxnID=str())
        # CASE 1: Ignore list
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="OXALODECARB-RXN",
            directory=dirBiocyc)
        pathways.createAndCheckSinksForRxn(
            model=testModel, rxnID="OXALODECARB_RXN_c",
            ignoreList=["PROTON_c", "PYRUVATE_c"])
        toCheck = [
            "CARBON_DIOXIDE_c", "OXALACETIC_ACID_c"]
        self.assertTrue(all([
            testSink in [
                sink.id for sink in testModel.sinks] for testSink in [
                    f'SK_{testMeta}' for testMeta in toCheck]
        ]))
        # CASE 2: no ignoreList
        toCheck = ["PYRUVATE_c", "OXALACETIC_ACID_c"]
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="OXALODECARB-RXN",
            directory=dirBiocyc)
        pathways.createAndCheckSinksForRxn(
            model=testModel, rxnID="OXALODECARB_RXN_c")
        self.assertTrue(all([
            testSink in [
                sink.id for sink in testModel.sinks] for testSink in [
                    f'SK_{testMeta}' for testMeta in toCheck]
        ]))
    @unittest.skip
    def test_testReactionForSolution(self):
        # CASE 1: Single Regular reaction
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="RXN-2206", directory=dirBiocyc)
        pathways.testReactionForSolution(
            model=testModel, rxnID="RXN_2206_c")
        self.assertGreater(abs(testModel.slim_optimize()), 0)
        # CASE 2: direction right to left, with ignoreList
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="1.8.4.9-RXN", directory=dirBiocyc)
        pathways.testReactionForSolution(
            model=testModel, rxnID="1.8.4.9_RXN_c",
            solutionRange=(0.01, 1000))
        self.assertGreater(abs(testModel.slim_optimize()), 0)
        # CASE 3: single reaction with ignoreList. i.e two reactions, in which
        # one metabolite X is created separately, and will be ignore in the
        # new function
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="RXN-2206", directory=dirBiocyc)
        testModel.add_boundary(
            testModel.metabolites.get_by_id("OXYGEN_MOLECULE_c"), "sink")
        pathways.testReactionForSolution(
            model=testModel, rxnID="RXN_2206_c",
            ignoreList=["OXYGEN_MOLECULE_c"])
        self.assertGreater(abs(testModel.slim_optimize()), 0)
        # CASE 4: single reverse reaction (as CASE 2) with a ignoreList
        testModel = cb.Model(0)
        add_reaction_from_root(
            model=testModel, xmlRoot="1.8.4.9-RXN", directory=dirBiocyc)
        testModel.add_boundary(
            testModel.metabolites.get_by_id("PROTON_c"), "sink")
        pathways.testReactionForSolution(
            model=testModel, rxnID="1.8.4.9_RXN_c",
            ignoreList=["PROTON_c"])
        self.assertGreater(abs(testModel.slim_optimize()), 0)
        # CASE 5: Stacking reactions: belongs to test_addPathModel
    @unittest.skip
    def test_addPathtoModel(self):
        # CASE 0a: wrong Model
        self.assertRaises(
            TypeError, pathways.addPathtoModel, model="Model",
            listReactions=list())
        # CASE 0b: wrong Reaction types
        self.assertRaises(
            TypeError, pathways.addPathtoModel, model=cb.Model(0),
            listReactions=[
                str(), tuple(), cb.Reaction(0)
            ])
        # CASE 1: Normal usage (3 reactions)
        testModel = cb.Model(0)
        add_meta_line_to_model(
            line="WATER, c", model=testModel, directory=dirBiocyc)
        add_meta_line_to_model(
            line="OXYGEN-MOLECULE, c", model=testModel, directory=dirBiocyc)
        for meta in ["WATER_c", "OXYGEN_MOLECULE_c"]:
            testModel.add_boundary(
                testModel.metabolites.get_by_id(meta), "sink")
        add_reaction_line_to_model(
            line=(
                f'Redox_homoproteins_c, Redox_homoproteins_c |'
                f'Ox-NADPH-Hemoprotein-Reductases_c:-1, '
                f'Red-NADPH-Hemoprotein-Reductases_c: 1'),
            model=testModel, directory=dirBiocyc)
        testList = pathways.create_reactions_for_list(
            model=testModel,
            pathwayList=["RXN-2206", "RXN-11414", "RXN-11422"],
            directory=dirBiocyc)
        pathways.addPathtoModel(
            model=testModel, listReactions=testList, ignoreList=[
                "WATER_c", "OXYGEN_MOLECULE_c"])
        self.assertGreater(abs(testModel.slim_optimize()), 0)
        # CASE 2: Stopping when unbalanced (by default, some of the reactions
        # are unbalanced)
        testModel = cb.Model(0)
        add_meta_line_to_model(
            line="WATER, c", model=testModel, directory=dirBiocyc)
        add_meta_line_to_model(
            line="OXYGEN-MOLECULE, c", model=testModel, directory=dirBiocyc)
        for meta in ["WATER_c", "OXYGEN_MOLECULE_c"]:
            testModel.add_boundary(
                testModel.metabolites.get_by_id(meta), "sink")
        add_reaction_line_to_model(
            line=(
                f'Redox_homoproteins_c, Redox_homoproteins_c |'
                f'Ox-NADPH-Hemoprotein-Reductases_c:-1, '
                f'Red-NADPH-Hemoprotein-Reductases_c: 1'),
            model=testModel, directory=dirBiocyc)
        testList = pathways.create_reactions_for_list(
            model=testModel,
            pathwayList=["RXN-2206", "RXN-11414", "RXN-11422"],
            directory=dirBiocyc)
        self.assertRaises(
            Warning, pathways.addPathtoModel, model=testModel,
            listReactions=testList, ignoreList=["WATER_c"], stopIfWrongMB=True)
        # CASE 3: regular path
        testModel = cb.Model(0)
        add_meta_from_file(
            model=testModel, filename=dirInput.joinpath(
                "metaToAdd_04_pathway_test.txt"),
            directory=dirBiocyc)
        # These metabolites would have normally other reactions to consume,
        # but since this is is testing, these extra metabolites needs to be
        # ignored !! 
        addSink = [
            "PROTON_c", "WATER_c", "OXYGEN_MOLECULE_c",
            "3_5_ADP_c", "CO_A_c"]
        [testModel.add_boundary(
            testModel.metabolites.get_by_id(meta), "sink") for meta in addSink]
        testList = pathways.createPathway(
            xmlRoot=get_xml_from_biocyc(
                bioID="PWY-1187", directory=dirBiocyc))
        testRxns = pathways.create_reactions_for_list(
            pathwayList=testList[1], directory=dirBiocyc, model=testModel)
        add_reaction_line_to_model(
            line=(
                f'Redox_homoproteins_c, Redox_homoproteins_c |'
                f'Ox-NADPH-Hemoprotein-Reductases_c:-1, '
                f'Red-NADPH-Hemoprotein-Reductases_c: 1'),
            model=testModel, directory=dirBiocyc)
        pathways.addPathtoModel(
            model=testModel, listReactions=testRxns,
            ignoreList=addSink)
        self.assertGreater(abs(testModel.slim_optimize()), 0)
        # CASE 4: Stacking 2 reactions
        testRxns = pathways.create_reactions_for_list(
            pathwayList=testList[2], directory=dirBiocyc, model=testModel)
        pathways.addPathtoModel(
            model=testModel, listReactions=testRxns,
            ignoreList=addSink)
        self.assertGreater(abs(testModel.slim_optimize()), 0)
    @unittest.skip
    def test_testAndAddCompletePathway(self):
        # CASE 1:
        testModel = cb.Model(0)
        add_meta_from_file(
            model=testModel, filename=dirInput.joinpath(
                "metaToAdd_04_pathway_test.txt"),
            directory=dirBiocyc)
        addSink = [
            "PROTON_c", "WATER_c", "OXYGEN_MOLECULE_c",
            "3_5_ADP_c", "CO_A_c", "CARBON_DIOXIDE_c", "CPD_3746_c"]
        [testModel.add_boundary(
            testModel.metabolites.get_by_id(meta), "sink") for meta in addSink]
        add_reaction_line_to_model(
            line=(
                f'Redox_homoproteins_c, Redox_homoproteins_c |'
                f'Ox-NADPH-Hemoprotein-Reductases_c:-1, '
                f'Red-NADPH-Hemoprotein-Reductases_c: 1'),
            model=testModel, directory=dirBiocyc)
        pathways.testAndAddCompletePathway(
            model=testModel, xmlRoot="PWY-1187", directory=dirBiocyc,
            ignoreList=addSink
        )
        # CASE 2: stacking another pathways (independent from each other)
        pathways.testAndAddCompletePathway(
            model=testModel, xmlRoot="AMMOXID-PWY", directory=dirBiocyc,
            ignoreList=addSink
        )
        self.assertGreater(abs(testModel.slim_optimize()), 0)
        sol = testModel.optimize()
        for demand in testModel.demands:
            self.assertGreater(abs(sol.fluxes[demand.id]), 0)
        testModel = cb.Model(0)
        add_meta_from_file(
            model=testModel, filename=dirInput.joinpath(
                "metaToAdd_05_pathway_test.txt"),
            directory=dirBiocyc)
        addSink = [
            "PROTON_c", "WATER_c", "OXYGEN_MOLECULE_c", "CO_A_c",
            "CARBON_DIOXIDE_c", "MAL_c", "NAD_c", "NADH_c"]
        [testModel.add_boundary(
            testModel.metabolites.get_by_id(meta), "sink") for meta in addSink]
        pathways.testAndAddCompletePathway(
            model=testModel, xmlRoot="PWY-5690", directory=dirBiocyc,
            ignoreList = addSink)
        # NOTE: check for correct demands


if __name__ == "__main__":
    unittest.main(verbosity=2)
