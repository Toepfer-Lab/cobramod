#!/usr/bin/env python3
from pathlib import Path
import unittest
import cobra as cb
from cobramod import pathways as pt
from cobramod.creation import get_xml_from_biocyc, add_reaction_from_root,\
    add_meta_from_string, add_meta_from_file, _add_reaction_line_to_model
from itertools import chain

dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_biocyc = Path.cwd().joinpath("tests").joinpath("data").joinpath("biocyc")

if not dir_biocyc.exists():
    dir_biocyc.mkdir(parents=True)


class ModulTesting(unittest.TestCase):

    def test__create_vertex_dict(self):
        self.assertRaises(
            TypeError, pt._create_vertex_dict, root=str())
        test_dict, _ = pt._create_vertex_dict(
                get_xml_from_biocyc(
                    identifier="PWY-1187", directory=dir_biocyc,
                    database="META"))
        self.assertIsInstance(test_dict, dict)
        # TODO: check for correct dictionary

    def test__find_end_vertex(self):
        # CASE 1: regular lineal / 2 child edges
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E"}
        self.assertIn("E", pt._find_end_vertex(vertex_dict=test_dict))
        # CASE 2: lineal with a cycle in the middle
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G"}
        self.assertIn("G", pt._find_end_vertex(vertex_dict=test_dict))
        # CASE 3: cyclical (regular), other cases will be examined further
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A"
            }
        self.assertRaises(
            StopIteration, next, pt._find_end_vertex(
                vertex_dict=test_dict))

    def test__find_start_vertex(self):
        """" Start vertex are always keys """
        # CASE 1: regular lineal
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E"}
        self.assertIn("A", pt._find_start_vertex(vertex_dict=test_dict))
        # CASE 2: lineal with 2 paths
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G"}
        self.assertIn("A", pt._find_start_vertex(vertex_dict=test_dict))
        self.assertIn("D", pt._find_start_vertex(vertex_dict=test_dict))
        # CASE 3: cyclical (regular), other cases will be examined further
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A"
            }
        self.assertRaises(
            StopIteration, next, pt._find_start_vertex(
                vertex_dict=test_dict))

    def test__verify_return_end_vertex(self):
        # CASE 1: regular lineal
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E"}
        self.assertIn(
            "E", pt._verify_return_end_vertex(vertex_dict=test_dict))
        # CASE 2: lineal with a cycle in the middle
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G"}
        self.assertIn(
            "G", pt._verify_return_end_vertex(vertex_dict=test_dict))
        # CASE 3: cyclical (regular), other cases will be examined further
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A"
            }
        self.assertEqual(
            1, len(
                list(pt._verify_return_end_vertex(vertex_dict=test_dict)))
        )
        # CASE 4a: cyclical with reactions in between
        test_dict = {
            # cycle
            "A": "B",
            "B": "C",  # could also be G
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F"
        }
        self.assertEqual(
            1, len(
                list(pt._verify_return_end_vertex(vertex_dict=test_dict)))
        )
        # CASE 4b: cyclical with reactions in between (G instead of C)
        test_dict = {
            # cycle
            "A": "B",
            "B": "G",  # could also be C
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F"
        }
        self.assertEqual(
            1, len(
                list(pt._verify_return_end_vertex(vertex_dict=test_dict)))
        )
        # CASE 5: cyclycal with one extra out
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "E",
            "D": "B",
            "E": "D"
            }
        self.assertEqual(
            1, len(
                list(pt._verify_return_end_vertex(vertex_dict=test_dict)))
        )

    def test_get_graph(self):
        # CASE 1: regular lineal
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E"}
        test_graph = pt.get_graph(vertex_dict=test_dict)
        self.assertIn(
            # Real path direction ->
            ["E", "D", "C", "B", "A"], [
                list(sequence) for sequence in test_graph])
        # CASE 2: lineal with a cycle in the middle
        test_dict = {
            "A": "B",
            "B": "C",  # D
            "C": "F",
            "D": "E",
            "E": "F",
            "F": "G"}
        test_graph = pt.get_graph(vertex_dict=test_dict)
        # 2 sequences
        self.assertEqual(2, len(list(test_graph)))
        # CASE 3: cyclical (reg,lar)
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A"
            }
        test_graph = pt.get_graph(vertex_dict=test_dict)
        self.assertEqual(
            6, len([x for x in set(chain.from_iterable(test_graph))])
            )
        # CASE 4a: cyclical with reactions in between
        # This can only be fix with the next step, checking that all edges
        # are mentioned at least on in the graph

        # CASE 4b: as 4a but replacing value of B
        test_dict = {
            # cycle
            "A": "B",
            "B": "G",
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F"
        }
        test_graph = pt.get_graph(vertex_dict=test_dict)
        self.assertEqual(
            8, len([x for x in set(chain.from_iterable(test_graph))])
            )
        # CASE 5: cyclycal with one extra out
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "E",
            "D": "B",
            "E": "D"
            }
        test_graph = pt.get_graph(vertex_dict=test_dict)
        self.assertEqual(
            5, len([x for x in set(chain.from_iterable(test_graph))])
            )

    def test__return_verified_graph(self):
        # CASE 0: regular lineal graph
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E"}
        test_graph = list(
            pt._return_verified_graph(
                vertex_dict=test_dict, vertex_set=set(
                    chain.from_iterable(test_dict.items())
                    )))
        self.assertIn(
            ["E", "D", "C", "B", "A"], [
                list(sequence) for sequence in test_graph])
        # CASE 1: cyclical (regular)
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E",
            "E": "F",
            "F": "A"
            }
        test_graph = list(
            pt._return_verified_graph(
                vertex_dict=test_dict, vertex_set=set(
                    chain.from_iterable(test_dict.items())
                    )))
        self.assertEqual(
            6, len([x for x in set(chain.from_iterable(test_graph))])
            )
        # CASE 3: cyclical with reactions in between
        test_dict = {
            # cycle
            "A": "B",
            "B": "C",  # could also be G
            "C": "Z",
            "Z": "D",
            "D": "E",
            "E": "A",
            # extra reactions
            "F": "Z",
            "G": "F"
        }
        test_graph = list(
            pt._return_verified_graph(
                vertex_dict=test_dict, vertex_set=set(
                    chain.from_iterable(test_dict.items())
                    )))
        # CASE 4: using replacement dict. As Key and value
        test_dict = {
            "A": "B",
            "B": "C",
            "C": "D",
            "D": "E"}
        test_replacement = {"A": "Z", "E": "F"}
        test_graph = list(
            pt._return_verified_graph(
                vertex_dict=test_dict,
                vertex_set=set(chain.from_iterable(test_dict.items())),
                replacement_dict=test_replacement
                    ))
        self.assertIn("Z", set(chain.from_iterable(test_graph)))
        self.assertIn("F", set(chain.from_iterable(test_graph)))

    def test_return_graph_from_root(self):
        # CASE 0: root not string or ET.Element
        self.assertRaises(
            AttributeError, pt.return_graph_from_root, root=int)
        # CASE 1: regular BioCycID
        test_list = pt.return_graph_from_root(
            root="PWY-1187", directory=dir_biocyc, database="META")
        # 3 end-metabolites + 1 in-between path (cyclic)
        self.assertEqual(len(
            set(chain.from_iterable(test_list))), 14)
        # CASE 2: cycle (no repetitions) / root is not a str
        test_list = pt.return_graph_from_root(
            root=get_xml_from_biocyc(
                identifier="PWY-5690", directory=dir_biocyc,
                database="META"))
        # Only one
        self.assertEqual(len(
            set(chain.from_iterable(test_list))), 9)
        # CASE 2: Cyclical with edges inbetween
        test_list = pt.return_graph_from_root(
            root="CALVIN-PWY", directory=dir_biocyc, database="META")
        # Total of 12, withouth repetitions
        self.assertEqual(len(
            set(chain.from_iterable(test_list))), 13)

    def test__create_reactions_for_iter(self):
        # CASE 1: Regular pathway
        test_model = cb.Model(0)
        test_graph = list(
            pt.return_graph_from_root(
                # its only 1 path
                root="PWY-1187", directory=dir_biocyc, database="META"))
        # Objects
        for sequence in test_graph:
            test_sequence_objects = pt._create_reactions_for_iter(
                sequence=sequence, directory=dir_biocyc, model=test_model,
                database="META")
            self.assertTrue(all(
                (isinstance(rxn, cb.Reaction) for rxn in test_sequence_objects)
                ))
        # Total amount
        self.assertEqual(14, len(
            set(chain.from_iterable(test_graph))
        ))
        # CASE 2: cycle (no repetitions)
        test_model = cb.Model(0)
        test_graph = list(pt.return_graph_from_root(
              # only 1 path
            root="PWY-5690", directory=dir_biocyc, database="META"))
        # Objects
        for sequence in test_graph:
            test_sequence_objects = pt._create_reactions_for_iter(
                sequence=sequence, directory=dir_biocyc, model=test_model,
                database="META")
            self.assertTrue(all(
                (isinstance(rxn, cb.Reaction) for rxn in test_sequence_objects)
                ))
        # Total amount
        self.assertEqual(9, len(
            set(chain.from_iterable(test_graph))
        ))
        # CASE 3: cycle with edgeds inbetween
        test_model = cb.Model(0)
        test_graph = list(pt.return_graph_from_root(
            root="CALVIN-PWY", directory=dir_biocyc,
            database="META"))  # only 1 path
        # Objects
        for sequence in test_graph:
            test_sequence_objects = pt._create_reactions_for_iter(
                sequence=sequence, directory=dir_biocyc, model=test_model,
                database="META")
            self.assertTrue(all(
                (isinstance(rxn, cb.Reaction) for rxn in test_sequence_objects)
                ))
        # Total amount
        self.assertEqual(13, len(
            set(chain.from_iterable(test_graph))
        ))
        # CASE 4: Compartment p
        test_model = cb.Model(0)
        test_model.compartments = {
            "e": "extracellular",
            "c": "cytosol",
            "p": "plastid"}
        test_graph = list(pt.return_graph_from_root(
            root="PWY-1187", directory=dir_biocyc, database="META"))
        for sequence in test_graph:
            test_sequence_objects = pt._create_reactions_for_iter(
                sequence=sequence, directory=dir_biocyc, model=test_model,
                compartment="p", database="META")
            (self.assertTrue(
                "p" in rxn.compartments) for rxn in test_sequence_objects)
        # CASE 5: same pathways with different compartment
        test_graph = list(pt.return_graph_from_root(
            root="PWY-1187", directory=dir_biocyc, database="META"))
        for sequence in test_graph:
            test_sequence_objects = pt._create_reactions_for_iter(
                sequence=sequence, directory=dir_biocyc, model=test_model,
                compartment="c", database="META")
            (self.assertTrue(
                "c" in rxn.compartments) for rxn in test_sequence_objects)

    def test__find_next_demand(self):
        # CASE 0a: invalid Model
        self.assertRaises(
            TypeError, pt._find_next_demand, model="Model",
            reaction_id=str())
        # CASe 0b: no metabolites found
        test_model = cb.Model(0)
        test_reaction = cb.Reaction(id="test")
        test_model.add_reactions([test_reaction])
        self.assertRaises(
            Warning, pt._find_next_demand, model=test_model,
            reaction_id="test")
        # creating reaction and model
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="OXALODECARB-RXN",
            directory=dir_biocyc, database="META")
        # test_model.add_reactions([test_reaction])
        # CASE 1a: left --> right
        self.assertEqual(
            "PYRUVATE_c",
            pt._find_next_demand(
                model=test_model, reaction_id="OXALODECARB_RXN_c")
        )
        # CASE 1b: with a ignore_list
        self.assertEqual(
            "CARBON_DIOXIDE_c",
            pt._find_next_demand(
                model=test_model, reaction_id="OXALODECARB_RXN_c",
                ignore_list=["PYRUVATE_c"])
        )
        # CASE 2a: left <-- right
        test_model.reactions.get_by_id("OXALODECARB_RXN_c").bounds = (-1000, 0)
        self.assertEqual(
            "OXALACETIC_ACID_c",
            pt._find_next_demand(
                model=test_model, reaction_id="OXALODECARB_RXN_c")
        )
        # CASE 2b: with a ignore_list
        self.assertEqual(
            "PROTON_c",
            pt._find_next_demand(
                model=test_model, reaction_id="OXALODECARB_RXN_c",
                ignore_list=["OXALACETIC_ACID_c", "FAKEID"])
        )
        # NOTE: Test for reversible reactions

    def test__fix_meta_for_boundaries(self):
        # CASE 0: ignore List
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="OXALODECARB-RXN",
            directory=dir_biocyc, database="META")
        self.assertRaises(
            Warning, pt._fix_meta_for_boundaries, model=test_model,
            metabolite="PROTON_c", ignore_list=["PROTON_c"])
        # CASE 1: PROTON, two normal reactions, one after another
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="OXALODECARB-RXN",
            directory=dir_biocyc, database="META")
        pt._fix_meta_for_boundaries(
            model=test_model, metabolite="PROTON_c")
        # Second reactions.
        add_reaction_from_root(
            model=test_model, root="AROMATIC-L-AMINO-ACID-DECARBOXYLASE-RXN",
            directory=dir_biocyc, database="META")
        pt._fix_meta_for_boundaries(
            model=test_model, metabolite="PROTON_c")
        self.assertNotIn(
            "SK_PROTON_c", [
                reaction.id for reaction in test_model.metabolites.get_by_id(
                    "PROTON_c").reactions])
        # CASE 2: Normal reaction, plus demand, plus test for sink
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="OXALODECARB-RXN",
            directory=dir_biocyc, database="META")
        test_model.add_boundary(
            test_model.metabolites.get_by_id("PROTON_c"), "demand")
        pt._fix_meta_for_boundaries(
            model=test_model, metabolite="PROTON_c")
        self.assertNotIn(
            "DM_PROTON_c", [
                reaction.id for reaction in test_model.metabolites.get_by_id(
                    "PROTON_c").reactions])
        # CASE 3: Adding an extra reaction, nothing should be left but
        # reactions
        add_reaction_from_root(
            model=test_model, root="AROMATIC-L-AMINO-ACID-DECARBOXYLASE-RXN",
            directory=dir_biocyc, database="META")
        pt._fix_meta_for_boundaries(
            model=test_model, metabolite="PROTON_c")
        self.assertNotIn(
            "SK_PROTON_c", [
                reaction.id for reaction in test_model.metabolites.get_by_id(
                    "PROTON_c").reactions])
        self.assertNotIn(
            "DM_PROTON_c", [
                reaction.id for reaction in test_model.metabolites.get_by_id(
                    "PROTON_c").reactions])

    def test__verify_side_sinks_for_rxn(self):
        # CASE 0: invalid Model
        self.assertRaises(
            TypeError, pt._verify_side_sinks_for_rxn, model=str(),
            rxn_id=str())
        self.assertRaises(
            ValueError, pt._verify_side_sinks_for_rxn,
            model=cb.Model(0),
            rxn_id=str(), side="Not")
        # CASE 1: normal creation (left side)
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="OXALODECARB-RXN", directory=dir_biocyc,
            database="META")
        pt._verify_side_sinks_for_rxn(
            model=test_model, rxn_id="OXALODECARB_RXN_c", side="left",
            ignore_list=["OXALACETIC_ACID_c"])
        self.assertTrue(
            "SK_PROTON_c" in (sink.id for sink in test_model.sinks)
        )
        # # CASE 2: Already sink
        pt._verify_side_sinks_for_rxn(
            model=test_model, rxn_id="OXALODECARB_RXN_c", side="left",
            ignore_list=["OXALACETIC_ACID_c"])
        # # CASE 3: normal creation (right side)
        pt._verify_side_sinks_for_rxn(
            model=test_model, rxn_id="OXALODECARB_RXN_c", side="right")
        self.assertTrue(
            "SK_PYRUVATE_c" in (sink.id for sink in test_model.sinks)
        )

    def test_verify_sinks_for_rxn(self):
        # CASE 0: Model invalid
        self.assertRaises(
            TypeError, pt.verify_sinks_for_rxn, model=str(),
            rxn_id=str())
        # CASE 1: Ignore list
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="OXALODECARB-RXN",
            directory=dir_biocyc, database="META")
        pt.verify_sinks_for_rxn(
            model=test_model, rxn_id="OXALODECARB_RXN_c",
            ignore_list=["PROTON_c", "PYRUVATE_c"])
        toCheck = [
            "CARBON_DIOXIDE_c", "OXALACETIC_ACID_c"]
        self.assertTrue(all([
            testSink in [
                sink.id for sink in test_model.sinks] for testSink in [
                    f'SK_{testMeta}' for testMeta in toCheck]
        ]))
        # CASE 2: no ignore_list
        toCheck = ["PYRUVATE_c", "OXALACETIC_ACID_c"]
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="OXALODECARB-RXN",
            directory=dir_biocyc, database="META")
        pt.verify_sinks_for_rxn(
            model=test_model, rxn_id="OXALODECARB_RXN_c")
        self.assertTrue(all([
            testSink in [
                sink.id for sink in test_model.sinks] for testSink in [
                    f'SK_{testMeta}' for testMeta in toCheck]
        ]))

    def test__test_rxn_for_solution(self):
        # CASE 1: Single Regular reaction
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="RXN-2206", directory=dir_biocyc,
            database="META")
        test_model.objective = "RXN_2206_c"
        pt._test_rxn_for_solution(
            model=test_model, rxn_id="RXN_2206_c")
        self.assertGreater(abs(test_model.slim_optimize()), 0)
        # CASE 2: direction right to left, with ignore_list
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="1.8.4.9-RXN", directory=dir_biocyc,
            database="META")
        test_model.objective = "1.8.4.9_RXN_c"
        test_model.objective_direction = "min"
        pt._test_rxn_for_solution(
            model=test_model, rxn_id="1.8.4.9_RXN_c",
            solutionRange=(0.01, 1000))
        self.assertGreater(abs(test_model.slim_optimize()), 0)
        # CASE 3: single reaction with ignore_list. i.e two reactions, in which
        # one metabolite X is created separately, and will be ignore in the
        # new function
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="RXN-2206", directory=dir_biocyc,
            database="META")
        test_model.add_boundary(
            test_model.metabolites.get_by_id("OXYGEN_MOLECULE_c"), "sink")
        test_model.objective = "RXN_2206_c"
        pt._test_rxn_for_solution(
            model=test_model, rxn_id="RXN_2206_c",
            ignore_list=["OXYGEN_MOLECULE_c"])
        self.assertGreater(abs(test_model.slim_optimize()), 0)
        # CASE 4: single reverse reaction (as CASE 2) with a ignore_list
        test_model = cb.Model(0)
        add_reaction_from_root(
            model=test_model, root="1.8.4.9-RXN", directory=dir_biocyc,
            database="META")
        test_model.add_boundary(
            test_model.metabolites.get_by_id("PROTON_c"), "sink")
        test_model.objective = "1.8.4.9_RXN_c"
        test_model.objective_direction = "min"
        pt._test_rxn_for_solution(
            model=test_model, rxn_id="1.8.4.9_RXN_c",
            ignore_list=["PROTON_c"])
        self.assertGreater(abs(test_model.slim_optimize()), 0)
        # CASE 5: Stacking reactions: belongs to test_addPathModel

    def test__add_sequence(self):
        # CASE 0a: wrong Model
        self.assertRaises(
            TypeError, pt._add_sequence, model="Model",
            sequence=list())
        # CASE 0b: wrong Reaction types
        self.assertRaises(
            TypeError, pt._add_sequence, model=cb.Model(0),
            sequence=[
                str(), tuple(), cb.Reaction(0)
            ])
        # CASE 1: Normal usage (3 reactions)
        test_model = cb.Model(0)
        add_meta_from_string(
            line="WATER, c", model=test_model, directory=dir_biocyc,
            database="META")
        add_meta_from_string(
            line="OXYGEN-MOLECULE, c", model=test_model, directory=dir_biocyc,
            database="META")
        for meta in ["WATER_c", "OXYGEN_MOLECULE_c"]:
            test_model.add_boundary(
                test_model.metabolites.get_by_id(meta), "sink")
        _add_reaction_line_to_model(
            line=(
                "Redox_homoproteins_c, Redox_homoproteins_c |"
                "Ox-NADPH-Hemoprotein-Reductases_c:-1,"
                "Red-NADPH-Hemoprotein-Reductases_c: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        # Dummy objective reaction
        _add_reaction_line_to_model(
            line=(
                "Dummy_c, Dummy_c reaction |" +
                "WATER_c: -1, OXYGEN_MOLECULE_c: 1"),
            model=test_model, directory=dir_biocyc)
        test_model.objective = "Dummy_c"
        test_list = list(pt._create_reactions_for_iter(
            model=test_model,
            sequence=["RXN-2206", "RXN-11414", "RXN-11422"],
            directory=dir_biocyc, database="META"))
        pt._add_sequence(
            model=test_model, sequence=test_list, ignore_list=[
                "WATER_c", "OXYGEN_MOLECULE_c"])
        self.assertGreater(abs(test_model.slim_optimize()), 0)
        # CASE 2: Stopping when unbalanced (by default, some of the reactions
        # are unbalanced)
        test_model = cb.Model(0)
        add_meta_from_string(
            line="WATER, c", model=test_model, directory=dir_biocyc,
            database="META")
        add_meta_from_string(
            line="OXYGEN-MOLECULE, c", model=test_model, directory=dir_biocyc,
            database="META")
        for meta in ["WATER_c", "OXYGEN_MOLECULE_c"]:
            test_model.add_boundary(
                test_model.metabolites.get_by_id(meta), "sink")
        # Dummy objective reaction
        _add_reaction_line_to_model(
            line=(
                "Dummy_c, Dummy_c reaction |" +
                "CARBON-DIOXIDE_c: -1, WATER_c: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        test_model.objective = "Dummy_c"
        _add_reaction_line_to_model(
            line=(
                "Redox_homoproteins_c, Redox_homoproteins_c |"
                "Ox-NADPH-Hemoprotein-Reductases_c:-1,"
                "Red-NADPH-Hemoprotein-Reductases_c: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        test_list = list(pt._create_reactions_for_iter(
            model=test_model,
            sequence=["RXN-2206", "RXN-11414", "RXN-11422"],
            directory=dir_biocyc, database="META"))
        # FIXME: reaction is not passing
        # self.assertRaises(
        #     Warning, pt._add_sequence, model=test_model,
        #     sequence=test_list, ignore_list=["WATER_c"], stop_wrong=True)
        # CASE 3: regular path
        test_model = cb.Model(0)
        add_meta_from_file(
            model=test_model, filename=dir_input.joinpath(
                "metaToAdd_04_pathway_test.txt"),
            directory=dir_biocyc, database="META")
        # These metabolites would have normally other reactions to consume,
        # but since this is is testing, these extra metabolites needs to be
        # ignored !
        new_sink_list = [
            "PROTON_c", "WATER_c", "OXYGEN_MOLECULE_c",
            "3_5_ADP_c", "CO_A_c"]
        for meta in new_sink_list:
            test_model.add_boundary(
                metabolite=test_model.metabolites.get_by_id(meta),
                type="sink")
        test_list = list(pt.return_graph_from_root(
            root="PWY-1187", directory=dir_biocyc, database="META"))
        test_reactions = list(
            pt._create_reactions_for_iter(
                sequence=test_list[1], directory=dir_biocyc,
                model=test_model, database="META"))
        _add_reaction_line_to_model(
            line=(
                "Redox_homoproteins_c, Redox_homoproteins_c |"
                "Ox-NADPH-Hemoprotein-Reductases_c:-1,"
                "Red-NADPH-Hemoprotein-Reductases_c: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        # Dummy objective reaction
        _add_reaction_line_to_model(
            line=(
                "Dummy_c, Dummy_c reaction |" +
                "CARBON-DIOXIDE_c: -1, WATER_c: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        test_model.objective = "Dummy_c"
        pt._add_sequence(
            model=test_model, sequence=test_reactions,
            ignore_list=new_sink_list, database="META")
        self.assertGreater(abs(test_model.slim_optimize()), 0)
        # CASE 4: Stacking 2 reactions
        test_reactions = list(
            pt._create_reactions_for_iter(
                sequence=test_list[2], directory=dir_biocyc, model=test_model,
                database="META"))
        pt._add_sequence(
            model=test_model, sequence=test_reactions,
            ignore_list=new_sink_list)
        self.assertGreater(abs(test_model.slim_optimize()), 0)

    def test__add_graph_from_root(self):
        # CASE 1:
        test_model = cb.Model(0)
        add_meta_from_file(
            model=test_model, filename=dir_input.joinpath(
                "metaToAdd_04_pathway_test.txt"),
            directory=dir_biocyc, database="META")
        new_sink_list = [
            "PROTON_c", "WATER_c", "OXYGEN_MOLECULE_c",
            "3_5_ADP_c", "CO_A_c", "CARBON_DIOXIDE_c", "CPD_3746_c"]
        for meta in new_sink_list:
            test_model.add_boundary(
                metabolite=test_model.metabolites.get_by_id(meta),
                type="sink")
        _add_reaction_line_to_model(
            line=(
                "Redox_homoproteins_c, Redox_homoproteins_c |"
                "Ox-NADPH-Hemoprotein-Reductases_c:-1,"
                "Red-NADPH-Hemoprotein-Reductases_c: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        # Dummy objective reaction
        _add_reaction_line_to_model(
            line=(
                "Dummy_c, Dummy reaction |" +
                "CARBON_DIOXIDE_c: -1, OXYGEN_MOLECULE_c: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        test_model.objective = "Dummy_c"
        pt._add_graph_from_root(
            model=test_model, root="PWY-1187", directory=dir_biocyc,
            ignore_list=new_sink_list, database="META"
        )
        # CASE 2: stacking another pathways (independent from each other)
        pt._add_graph_from_root(
            model=test_model, root="AMMOXID-PWY", directory=dir_biocyc,
            ignore_list=new_sink_list, database="META"
        )
        self.assertGreater(abs(test_model.slim_optimize()), 0)
        sol = test_model.optimize()
        for demand in test_model.demands:
            self.assertGreater(abs(sol.fluxes[demand.id]), 0)
        # CASE 3: Same with case 1 but with another compartment
        test_model = cb.Model(0)
        test_model.compartments = {
            "e": "extracellular",
            "c": "cytosol",
            "p": "plastid"}
        add_meta_from_file(
            model=test_model, filename=dir_input.joinpath(
                "metaToAdd_04_pathway_p_comp_test.txt"),
            directory=dir_biocyc, database="META")
        new_sink_list = [
            "PROTON_p", "WATER_p", "OXYGEN_MOLECULE_p",
            "3_5_ADP_p", "CO_A_p", "CARBON_DIOXIDE_p", "CPD_3746_p"]
        for meta in new_sink_list:
            test_model.add_boundary(
                metabolite=test_model.metabolites.get_by_id(meta),
                type="sink")
        _add_reaction_line_to_model(
            line=(
                "Redox_homoproteins_c, Redox_homoproteins_p |"
                "Ox-NADPH-Hemoprotein-Reductases_p:-1,"
                "Red-NADPH-Hemoprotein-Reductases_p: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        # Dummy objective reaction
        _add_reaction_line_to_model(
            line=(
                "Dummy_c, Dummy reaction |" +
                "CARBON-DIOXIDE_p: -1, OXYGEN-MOLECULE_p: 1"),
            model=test_model, directory=dir_biocyc, database="META")
        test_model.objective = "Dummy_c"
        pt._add_graph_from_root(
            model=test_model, root="PWY-1187", directory=dir_biocyc,
            ignore_list=new_sink_list, compartment="p", database="META"
        )
        pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
