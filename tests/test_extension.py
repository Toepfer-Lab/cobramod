#!/usr/bin/env python3
"""Unittest for module extension

this module test the behavior of multiple functions, which are responsible to
add metabolites and reactions into a metabolic model. This module can be
divided in to parts:

- CreatingSequences: Functions, that create the corresponding reactions as
sequences and their corresponding flux test.
- AddingPathways: Functions, that manage the addition of Pathways into the
metabolic models.
"""

import logging
import unittest
from pathlib import Path

import cobra.core as cobra_core
import cobramod.error as cmod_error
import cobramod.retrieval as cmod_retrieval
import requests
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.core import extension as ex
from cobramod.core.creation import add_reactions
from cobramod.core.pathway import Pathway
from cobramod.debug import change_to_debug
from cobramod.parsing.db_version import DataVersionConfigurator
from cobramod.test import textbook, textbook_biocyc, textbook_kegg

NAME = "test_model"

change_to_debug()

dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class CreatingSequences(unittest.TestCase):
    @classmethod
    def setUp(cls):
        # ToDo Versioning should be tested in its own test class
        data_conf = DataVersionConfigurator()
        data_conf.ignore_db_versions = True

    def test_create_reactions(self):
        # CASE: Simple Case Biocyc
        generator = ex.yield_reaction_from_list(
            sequence=[
                "OXALODECARB-RXN",
                "AROMATIC-L-AMINO-ACID-DECARBOXYLASE-RXN",
            ],
            compartment="c",
            directory=dir_data,
            database="ARA",
            replacement={},
            show_imbalance=False,
            stop_imbalance=False,
            model=cobra_core.Model(NAME),
            model_id="",
            genome=None,
        )
        self.assertIsInstance(obj=next(generator), cls=cobra_core.Reaction)

        # CASE: Simple case Kegg
        generator = ex.yield_reaction_from_list(
            sequence=["R00200", "R00114"],
            compartment="p",
            directory=dir_data,
            database="KEGG",
            replacement={},
            show_imbalance=False,
            stop_imbalance=False,
            model=cobra_core.Model(NAME),
            model_id="",
            genome="eco",
        )
        self.assertIsInstance(obj=next(generator), cls=cobra_core.Reaction)
        # CASE: Showing when unbalanced
        generator = ex.yield_reaction_from_list(
            sequence=["RXN-11414", "RXN-11422"],
            compartment="c",
            directory=dir_data,
            database="ARA",
            replacement={},
            show_imbalance=True,
            stop_imbalance=False,
            model=cobra_core.Model(NAME),
            model_id="",
            genome=None,
        )
        # RXN-11414 is unbalanced
        with self.assertLogs(level=logging.DEBUG) as cm:
            next(generator)
            self.assertIn("unbalanced", cm.output[-1])

        # CASE: Stopping when unbalanced
        generator = ex.yield_reaction_from_list(
            sequence=["RXN-11414", "RXN-11422"],
            compartment="c",
            directory=dir_data,
            database="ARA",
            show_imbalance=False,
            stop_imbalance=True,
            replacement={},
            model=cobra_core.Model(NAME),
            model_id="",
            genome=None,
        )
        self.assertRaises(cmod_error.UnbalancedReaction, next, generator)

    def test_test_non_zero_flux(self):
        # CASE: Single Regular reaction
        test_model = cobra_core.Model(NAME)

        add_reactions(
            model=test_model,
            obj=["RXN-2206, c"],
            directory=dir_data,
            database="ARA",
            replacement={},
            show_imbalance=False,
        )
        ex.test_non_zero_flux(model=test_model, reaction="RXN_2206_c")
        self.assertEqual(first=7, second=len(test_model.sinks))

        # CASE: minimun range not reached
        # INFO: this extra steps are needed in order to make it fail
        test_model = textbook_biocyc.copy()
        add_reactions(
            model=test_model,
            obj=[
                "Redox_ADP_ATP_p, Redox_ADP_ATP_p | ADP_p <-> ATP_p",
                "TRANS_Pi_cp, Transport Phosphate_cp | Pi_c <-> Pi_p",
                "TRANS_GLUTATHIONE_cp, Transport GLUTATHIONE_cp |"
                + " GLUTATHIONE_c <-> GLUTATHIONE_p",
                "GLUTATHIONE-SYN-RXN, p",
            ],
            database="ARA",
            directory=dir_data,
            show_imbalance=False,
        )

        # CASE: direction right to left
        test_model = cobra_core.Model(NAME)
        add_reactions(
            model=test_model,
            obj=["1.8.4.9-RXN, c"],
            directory=dir_data,
            database="ARA",
        )
        ex.test_non_zero_flux(model=test_model, reaction="1.8.4.9_RXN_c")
        self.assertEqual(first=6, second=len(test_model.sinks))


class AddingPathways(unittest.TestCase):
    """
    Test for functions related to the addition of Pathways and their
    visualizations.
    """

    def test_add_reactions_to_Pathway(self):
        # CASE: Normal usage (3 reactions)
        test_model = textbook_biocyc.copy()
        test_list = list(
            ex.yield_reaction_from_list(
                sequence=["RXN-2206", "RXN-11414", "RXN-11422"],
                compartment="c",
                directory=dir_data,
                database="ARA",
                replacement={},
                show_imbalance=False,
                stop_imbalance=False,
                model=test_model,
                model_id="",
                genome=None,
            )
        )
        test_group = Pathway("test_group")

        ex.add_reactions_to_Pathway(
            model=test_model,
            pathway=test_group,
            sequence=test_list,
            ignore_list=["WATER_c", "OXYGEN_MOLECULE_c"],
        )
        self.assertGreater(
            abs(test_model.slim_optimize(error_value=0)),
            0,  # type: ignore
        )
        self.assertEqual(len(test_group.members), 3)

        # CASE: reactions already in model
        test_model = textbook.copy()
        reactions = [
            reaction
            for reaction in test_model.reactions
            if reaction.id in ("GAPD", "PGK", "PGM")
        ]
        test_group = Pathway("test_group")

        ex.add_reactions_to_Pathway(
            model=test_model,
            pathway=test_group,
            sequence=reactions,
            ignore_list=[],
        )
        self.assertEqual(len(test_group.members), 3)

        # CASE: KEGG
        test_model = textbook_kegg.copy()
        test_list = list(
            ex.yield_reaction_from_list(
                sequence=["R00894", "R00497"],
                compartment="c",
                directory=dir_data,
                database="KEGG",
                replacement={},
                show_imbalance=False,
                stop_imbalance=False,
                model=test_model,
                model_id="",
                genome="ath",
            )
        )
        test_group = Pathway("test_group")

        ex.add_reactions_to_Pathway(
            model=test_model,
            pathway=test_group,
            sequence=test_list,
            ignore_list=["WATER_c", "OXYGEN_MOLECULE_c"],
        )
        self.assertGreater(abs(test_model.slim_optimize(error_value=0)), 0)
        self.assertEqual(len(test_group.members), 2)

    def test_add_pathway_from_dict(self):
        data = cmod_retrieval.get_data(
            identifier="M00118",
            directory=dir_data,
            # debug_level=10,
            database="KEGG",
        )
        # CASE 1: regular test with KEGG
        test_model = textbook_kegg.copy()
        ex.add_pathway_from_data(
            model=test_model,
            group=None,
            data=data,
            directory=dir_data,
            database="KEGG",
            compartment="c",
            avoid_list=[],
            replacement={},
            ignore_list=[],
            show_imbalance=False,
            stop_imbalance=False,
            model_id="",
            genome="hsa",
        )
        self.assertIn(
            member="M00118",
            container=[group.id for group in test_model.groups],
        )
        for item in ["2729", "2937"]:
            self.assertIn(
                member=item, container=[gene.id for gene in test_model.genes]
            )
        # CASE 2: Using another label
        test_model = textbook_kegg.copy()
        ex.add_pathway_from_data(
            model=test_model,
            data=data,
            directory=dir_data,
            database="KEGG",
            group="ALTERNATIVE",
            compartment="c",
            avoid_list=[],
            replacement={},
            ignore_list=[],
            show_imbalance=False,
            stop_imbalance=False,
            model_id="",
            genome="hsa",
        )
        self.assertIn(
            member="ALTERNATIVE",
            container=[group.id for group in test_model.groups],
        )

    def test_add_pathway_from_list(self):
        # CASE: regular test
        test_model = textbook_biocyc.copy()
        ex.add_pathway_from_strings(
            model=test_model,
            identifier="test_group",
            compartment="c",
            sequence=["R01063", "R09084"],
            database="KEGG",
            directory=dir_data,
            avoid_list=[],
            replacement={},
            ignore_list=[],
            show_imbalance=False,
            stop_imbalance=False,
            model_id="",
            genome="mba",
        )
        self.assertTrue(test_model.groups.has_id("test_group"))

        for item in ["Mbar_A2189", "Mbar_A3564"]:
            self.assertTrue(test_model.genes.has_id(item))

    def test_add_pathway_from_file(self):
        test_model = cobra_core.Model(NAME)

        file = dir_input.joinpath("reactions_identifiers.txt")

        ex.add_pathway_from_file(
            model=test_model,
            file=file,
            replacement={"Not_GLC": "GLC"},
            identifier="test_group",
            show_imbalance=False,
            stop_imbalance=False,
            model_id="",
            database="META",
            directory=dir_data,
            ignore_list=[],
            genome=None,
        )
        for reaction in (
            "RXN_14462_c",
            "ACETALD_DEHYDROG_RXN_c",
            "AMONITRO_RXN_c",
        ):
            self.assertIn(
                member=reaction,
                container=[reaction.id for reaction in test_model.reactions],
            )

    def test_add_pathway(self):
        # CASE: Regular Biocyc
        test_model = textbook_biocyc.copy()

        ex.add_pathway(
            model=test_model,
            pathway="PWY-1187",
            compartment="c",
            directory=dir_data,
            database="ARA",
            ignore_list=[],
            show_imbalance=False,
        )
        self.assertIn(
            member="RXN_11438_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        for item in ["AT1G18590", "AT1G74090"]:
            self.assertIn(
                member=item, container=[gene.id for gene in test_model.genes]
            )
        # CASE: stacking another pathways (independent from each other)
        ex.add_pathway(
            model=test_model,
            pathway="AMMOXID-PWY",
            compartment="c",
            directory=dir_data,
            database="META",
            show_imbalance=False,
        )
        self.assertGreater(abs(test_model.slim_optimize(error_value=0)), 0)
        sol = test_model.optimize()

        # Test for demands, all should have a flux
        for demand in test_model.demands:
            value = sol.fluxes[demand.id]
            if not isinstance(value, float):
                raise TypeError("Value is not float")

            self.assertGreater(abs(value), 0)
        self.assertEqual(first=len(test_model.groups), second=2)

        # CASE: using a simple sequence
        test_model = textbook_biocyc.copy()
        test_sequence = ["RXN-2206", "RXN-11414", "RXN-11422", "RXN-11430"]
        ex.add_pathway(
            model=test_model,
            compartment="c",
            pathway=test_sequence,
            ignore_list=[],
            database="ARA",
            directory=dir_data,
            show_imbalance=False,
        )
        self.assertGreater(a=test_model.slim_optimize(), b=0)
        self.assertIn(
            member="RXN_11422_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE: Behavior if reaction was already in model
        test_model = textbook_biocyc.copy()

        # Adding reactions
        test_sequence = ["RXN-2206", "RXN-11414"]
        ex.add_pathway(
            model=test_model,
            group="old_reactions",
            compartment="c",
            pathway=test_sequence,
            ignore_list=[],
            database="META",
            directory=dir_data,
            show_imbalance=False,
        )
        test_sequence = ["RXN-2206", "RXN-11414", "RXN-11422", "RXN-11430"]
        ex.add_pathway(
            model=test_model,
            compartment="c",
            pathway=test_sequence,
            ignore_list=[],
            database="ARA",
            directory=dir_data,
            show_imbalance=False,
        )
        self.assertGreater(a=test_model.slim_optimize(), b=0)

        test_group = test_model.groups.get_by_id("custom_group")
        if not isinstance(test_group, Pathway):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        for identifier in ["RXN_2206_c", "RXN_11414_c", "RXN_11422_c"]:
            self.assertIn(
                member=identifier,
                container=[reaction.id for reaction in test_group.members],
            )
        # CASE: KEGG simple pathway
        test_model = textbook_kegg.copy()
        ex.add_pathway(
            model=test_model,
            pathway="M00118",
            database="KEGG",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
            genome="hsa",
        )
        self.assertGreater(a=test_model.slim_optimize(), b=0)
        self.assertIn(
            member="R00894_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE: KEGG lineal coplex pathway
        test_model = textbook_kegg.copy()
        ex.add_pathway(
            model=test_model,
            pathway="M00001",
            database="KEGG",
            directory=dir_data,
            compartment="c",
            ignore_list=[],
            avoid_list=[],
            show_imbalance=False,
            genome="hsa",
        )
        self.assertGreater(a=test_model.slim_optimize(), b=0)
        self.assertIn(
            member="R01063_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        self.assertIn(
            member="R01518_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # CASE: Check for translations in pathways
        test_model = textbook_kegg.copy()
        ex.add_pathway(
            model=test_model,
            pathway="SALVADEHYPOX-PWY",
            database="ECO",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
        )
        self.assertGreater(a=test_model.slim_optimize(), b=0)
        test_reaction = test_model.reactions.get_by_id("ADENODEAMIN_RXN_c")

        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        self.assertIn(
            member="C00080_c",
            container=[
                metabolite.id for metabolite in test_reaction.metabolites
            ],
        )
        # CASE: From sequence. Skipping already in Model. Stacking in one.
        test_model = textbook_kegg.copy()
        ex.add_pathway(
            model=test_model,
            pathway=["ACALD", "MALS"],
            database="BIGG",
            model_id="e_coli_core",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
            ignore_list=["R00228_c", "R00472_c"],
        )

        ex.add_pathway(
            model=test_model,
            pathway=["ADENODEAMIN-RXN"],
            database="ECO",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
        )
        test_reaction = test_model.reactions.get_by_id("ADENODEAMIN_RXN_c")
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Given object is not a valid COBRApy Reaction")

        test_group = test_model.groups.get_by_id("custom_group")
        if not isinstance(test_group, Pathway):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        self.assertEqual(first=len(test_group.members), second=3)
        self.assertGreater(a=test_model.slim_optimize(), b=0)
        self.assertNotIn(
            member="ACALD_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        self.assertNotIn(
            member="MALS_c",
            container=[reaction.id for reaction in test_model.reactions],
        )
        # WATER
        self.assertIn(
            member="C00001_c",
            container=[
                metabolite.id for metabolite in test_reaction.metabolites
            ],
        )
        # CASE: Testing avoid list
        test_model = textbook_kegg.copy()
        test_sequence = ["RXN-2206", "RXN-11414", "RXN-11422", "RXN-11430"]
        ex.add_pathway(
            model=test_model,
            pathway=test_sequence,
            database="ARA",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
            avoid_list=["RXN-2206"],
        )
        test_group = test_model.groups.get_by_id("custom_group")
        if not isinstance(test_group, Pathway):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        self.assertEqual(first=len(test_group.members), second=3)
        self.assertEqual(first=len(test_group.graph), second=3)

        # CASE: Stacking pathways in one pathways
        test_sequence = ["RXN-11438", "RXN-2208", "RXN-2209", "RXN-2221"]
        ex.add_pathway(
            model=test_model,
            pathway=test_sequence,
            database="ARA",
            directory=dir_data,
            compartment="c",
            show_imbalance=False,
        )
        test_group = test_model.groups.get_by_id("custom_group")
        if not isinstance(test_group, Pathway):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        self.assertEqual(first=len(test_group.members), second=7)
        self.assertEqual(first=len(test_group.graph), second=7)

        # CASE: Using replacement arguments to rename
        test_model = textbook_biocyc.copy()
        test_replacement = {
            "RXN-8052": "Reaction_8052",
            "RXN-2208": "Reaction_2208",
            "RXN-2223": "Reaction_2223",
        }
        ex.add_pathway(
            model=test_model,
            pathway="PWY-1187",
            compartment="c",
            directory=dir_data,
            database="ARA",
            replacement=test_replacement,
            show_imbalance=False,
        )
        test_group = test_model.groups.get_by_id("PWY-1187")
        if not isinstance(test_group, Pathway):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        for meta in test_replacement.values():
            identifier = f"{meta}_c"
            self.assertIn(
                member=identifier,
                container=[reaction.id for reaction in test_model.reactions],
            )
            self.assertIn(
                member=identifier,
                container=[reaction.id for reaction in test_group.members],
            )

        # CASE: Using replacement arguments to replace
        test_model = textbook_biocyc.copy()
        test_sequence = ["RXN-2206", "RXN-11414", "RXN-11422", "RXN-11430"]

        ex.add_pathway(
            model=test_model,
            compartment="c",
            pathway=test_sequence,
            ignore_list=[],
            database="ARA",
            directory=dir_data,
            replacement={"RXN-2206": "ACALDt"},
            show_imbalance=False,
        )

        test_group = test_model.groups.get_by_id("custom_group")
        if not isinstance(test_group, Pathway):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        self.assertIn(
            member="ACALDt_c",
            container=[reaction.id for reaction in test_group.members],
        )
        # CASE: Raising AttributeError
        test_model = textbook_biocyc.copy()

        self.assertRaises(
            requests.HTTPError,
            ex.add_pathway,
            model=test_model,
            compartment="c",
            pathway=dir_input.joinpath("reactions_identifiers.txt"),
            database="ARA",
            ignore_list=[],
            directory=dir_data,
            show_imbalance=False,
        )

        # CASE: Using Path
        test_model = textbook_biocyc.copy()

        ex.add_pathway(
            model=test_model,
            compartment="c",
            pathway=dir_input.joinpath("reactions_identifiers.txt"),
            database="META",
            ignore_list=[],
            directory=dir_data,
            show_imbalance=False,
        )

        test_group = test_model.groups.get_by_id("custom_group")
        if not isinstance(test_group, Pathway):
            raise TypeError("Given object is not a valid Cobramod Pathway")

        for reaction in (
            "RXN_14462_c",
            "ACETALD_DEHYDROG_RXN_c",
            "AMONITRO_RXN_c",
        ):
            self.assertIn(
                member=reaction,
                container=[reaction.id for reaction in test_group.members],
            )


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
