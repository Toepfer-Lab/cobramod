#!/usr/bin/env python3
"""Unit-test for model behavior

This module checks if a large and a small model behaves as intended. This
examples should simulate real cases. There are two test:

- ShortModel: This should utilize the textbook_biocyc model from cobramod.
- LargeModel: Uses a real GEM
"""

import unittest
from pathlib import Path

import cobra.core as cobra_core
import cobra.io as cobra_io
import cobramod.error as cmod_error
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.core import extension as ex
from cobramod.core.creation import add_reactions
from cobramod.debug import change_to_debug
from cobramod.parsing.db_version import DataVersionConfigurator
from cobramod.test import textbook_biocyc

change_to_debug()
dir_data = Path(__file__).resolve().parent.joinpath("data")

# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")

dir_input = Path(__file__).resolve().parent.joinpath("input")

# Large Model
main_model = cobra_io.read_sbml_model(
    str(dir_input.joinpath("test_model.sbml"))
)


class ShortModel(unittest.TestCase):
    @classmethod
    def setUp(cls):
        data_conf = DataVersionConfigurator()
        data_conf.ignore_db_versions = True

    def test_lineal_pathways(self):
        # For this test, Gluconeogenesis; L-aspartate and L-asparagine
        # biosynthesis, and Nicotine biosynthesis added to test for feasible
        # results.
        test_model = textbook_biocyc.copy()

        # Adding Gluconeogenesis
        ex.add_pathway(
            model=test_model,
            pathway="GLUCONEO-PWY",
            directory=dir_data,
            database="META",
            compartment="c",
            show_imbalance=False,
        )
        self.assertGreater(a=abs(test_model.optimize().objective_value), b=0)
        self.assertIn(
            member="GLUCONEO-PWY",
            container=[group.id for group in test_model.groups],
        )
        add_reactions(
            model=test_model,
            obj=[
                "Redox_Hemoprotein_c, Redox_Hemoprotein_c | "
                + "Ox-NADPH-Hemoprotein-Reductases_c <-> "
                + "Red-NADPH-Hemoprotein-Reductases_c"
            ],
            directory=dir_data,
            database="ARA",
        )

        # Adding Nicotine pathway
        ex.add_pathway(
            model=test_model,
            pathway="PWY-5316",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
            show_imbalance=False,
        )
        self.assertIn(
            member="PWY-5316",
            container=[group.id for group in test_model.groups],
        )

        # Adding aspartate and asparagine biosynthesis
        # Before 27.1 in Metacyc, this was a pathway
        self.assertRaises(
            cmod_error.SuperpathwayException,
            ex.add_pathway,
            model=test_model,
            pathway="ASPASN-PWY",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
            show_imbalance=False,
        )
        ex.add_pathway(
            model=test_model,
            pathway="ASPARAGINESYN-PWY",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
            show_imbalance=False,
            group="ASPASN-PWY",
        )
        ex.add_pathway(
            model=test_model,
            pathway="ASPARTATESYN-PWY",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
            show_imbalance=False,
            group="ASPASN-PWY",
        )
        self.assertGreater(a=abs(test_model.optimize().objective_value), b=0)
        self.assertIn(
            member="ASPASN-PWY",
            container=[group.id for group in test_model.groups],
        )

    def test_cyclical_pathways(self):
        test_model = textbook_biocyc.copy()

        # Adding Mannitol cycle
        ex.add_pathway(
            model=test_model,
            pathway="PWY-6531",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
            show_imbalance=False,
        )
        self.assertIn(
            member="PWY-6531",
            container=[group.id for group in test_model.groups],
        )

        # Adding glyoxylate cycle
        ex.add_pathway(
            model=test_model,
            pathway="GLYOXYLATE-BYPASS",
            directory=dir_data,
            database="ARA",
            compartment="c",
            ignore_list=[],
            show_imbalance=False,
        )
        self.assertGreater(a=abs(test_model.optimize().objective_value), b=0)
        self.assertIn(
            member="GLYOXYLATE-BYPASS",
            container=[group.id for group in test_model.groups],
        )

    def test_multi_compartment(self):
        test_model = textbook_biocyc.copy()
        # This files has transports between cytosol and plastids
        add_reactions(
            model=test_model,
            obj=dir_input.joinpath("test_multi_reactions.txt"),
            database="ARA",
            directory=dir_data,
        )

        ex.add_pathway(
            model=test_model,
            pathway="GLUTATHIONESYN-PWY",
            directory=dir_data,
            database="ARA",
            ignore_list=[],
            compartment="p",
            show_imbalance=False,
        )
        reaction = test_model.reactions.get_by_id("Biomass_Ecoli_core")
        if not isinstance(reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")

        metabolite = test_model.metabolites.get_by_id("GLUTATHIONE_p")
        if not isinstance(metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy metabolite")

        reaction.add_metabolites({metabolite: -1})

        self.assertIn(
            member="GLUTATHIONESYN-PWY",
            container=[group.id for group in test_model.groups],
        )

        self.assertGreater(abs(test_model.optimize().objective_value), 0)

        add_reactions(
            model=test_model,
            obj=["TRANS_GLY_cp, Transport GLY_cp | GLY_c <-> GLY_p"],
            directory=dir_data,
            database="ARA",
        )

        ex.add_pathway(
            model=test_model,
            pathway="PWY-1187",
            directory=dir_data,
            database="ARA",
            ignore_list=["PYRUVATE_c", "CO_A_c", "PROTON_c", "CPD_3746_c"],
            compartment="c",
            show_imbalance=False,
        )
        reaction = test_model.reactions.get_by_id("Biomass_Ecoli_core")
        if not isinstance(reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")

        metabolite_a = test_model.metabolites.get_by_id(
            "3_METHYLSULFINYLPROPYL_GLUCOSINOLATE_c"
        )
        if not isinstance(metabolite_a, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy metabolite")
        metabolite_b = test_model.metabolites.get_by_id("GLUTATHIONE_p")
        if not isinstance(metabolite_b, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy metabolite")

        reaction.add_metabolites({metabolite_a: -0.75, metabolite_b: 1})
        self.assertGreater(abs(test_model.optimize().objective_value), 0)

        # Lipid initalization
        ex.add_pathway(
            model=test_model,
            pathway="PWY-4381",
            directory=dir_data,
            database="ARA",
            ignore_list=["CO_A_p"],
            compartment="p",
            show_imbalance=False,
        )
        reaction = test_model.reactions.get_by_id("Biomass_Ecoli_core")
        if not isinstance(reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")

        metabolite = test_model.metabolites.get_by_id("Acetoacetyl_ACPs_p")
        if not isinstance(metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy metabolite")

        reaction.add_metabolites({metabolite: -1})

        self.assertGreater(test_model.slim_optimize(error_value=0), 0)


class LargeModel(unittest.TestCase):
    def test_lineal_pathways(self):
        """In this test, unrelated pathways are appended to the large model
        and tested
        """
        test_model = main_model.copy()

        # Autotroph environment
        test_model.reactions.Sucrose_tx.bounds = (0, 0)  # SUCROSE
        test_model.reactions.GLC_tx.bounds = (0, 0)  # GLUCOSE

        # # Photon uptake
        test_model.reactions.Photon_tx.bounds = (-1000, 250)

        # No ammonium
        test_model.reactions.NH4_tx.bounds = (0, 0)
        test_model.objective = "Biomass_tx"

        # Adding methionine metabolism
        ex.add_pathway(
            model=test_model,
            pathway="PWY-7614",
            directory=dir_data,
            database="META",
            compartment="p",
            show_imbalance=False,
        )

        # Optimizing methiin
        test_model.objective = "RXN_8908_p"

        self.assertGreater(test_model.optimize().fluxes["RXN_16201_p"], 0)
        self.assertGreater(test_model.slim_optimize(error_value=0), 0)

        # Adding stachyose biosynthesis
        test_model.objective = "Biomass_tx"

        ex.add_pathway(
            model=test_model,
            pathway="PWY-5337",
            directory=dir_data,
            database="ARA",
            compartment="c",
            show_imbalance=False,
        )

        # Optimizing raffinose
        test_model.objective = "2.4.1.82_RXN_c"

        self.assertGreater(test_model.optimize().fluxes["2.4.1.123_RXN_c"], 0)
        self.assertGreater(test_model.slim_optimize(error_value=0), 0)

        # Adding abscisic acid
        test_model.objective = "Biomass_tx"

        ex.add_pathway(
            model=test_model,
            pathway="PWY-695",
            directory=dir_data,
            database="ARA",
            compartment="p",
            show_imbalance=False,
        )

        # Optimizing abscisate
        test_model.objective = "1.2.3.14_RXN_p"

        self.assertGreater(test_model.optimize().fluxes["1.1.1.288_RXN_p"], 0)
        self.assertGreater(test_model.slim_optimize(error_value=0), 0)

        test_model.objective = "Biomass_tx"
        self.assertGreater(test_model.slim_optimize(error_value=0), 0)


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
