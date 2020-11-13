#!/usr/bin/env python3
from logging import DEBUG
from pathlib import Path
import unittest

import cobra as cb

from cobramod import extension as ex
from cobramod.creation import (
    _add_reaction_line_to_model,
    add_reactions_from_file,
)
from cobramod.debug import debug_log

debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")
# Short model
main_model1 = cb.io.read_sbml_model(
    str(dir_input.joinpath("test_model01.sbml"))
)
# Large Model
main_model2 = cb.io.read_sbml_model(
    str(dir_input.joinpath("test_model02.sbml"))
)


if not dir_data.exists():
    dir_data.mkdir(parents=True)


# TODO: add test with replacment_dicts
class TestingShortModel(unittest.TestCase):
    def test_appending_lineal_pathways(self):
        """For this test, Gluconeogenesis; L-aspartate and L-asparagine
        biosynthesis, and Nicotine biosynthesis added to test for feasible
        results. Due to a small model,  ATP needs to enter as a sink reaction
        but is later removed from model without interfering with the the new
        pathways
        """
        test_model = main_model1.copy()
        test_model.add_boundary(
            test_model.metabolites.get_by_id("ATP_c"), "sink"
        )
        # Adding Gluconeogenesis
        ex.add_graph_to_model(
            model=test_model,
            graph="GLUCONEO-PWY",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=["MAL_c", "PROTON_c", "Pi_c", "ATP_c"],
        )
        test_model.remove_reactions(["SK_ATP_c"])
        self.assertGreater(abs(test_model.optimize().objective_value), 0)
        _add_reaction_line_to_model(
            model=test_model,
            line=(
                "Redox_Hemoprotein_c, Redox_Hemoprotein_c | "
                + "Ox-NADPH-Hemoprotein-Reductases_c:-1, "
                + "Red-NADPH-Hemoprotein-Reductases_c: 1"
            ),
            directory=dir_data,
            database="META",
        )
        # Adding Nicotine pathway
        ex.add_graph_to_model(
            model=test_model,
            graph="PWY-5316",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
        )
        # Ading aspartate and asparagine biosynthesis
        ex.add_graph_to_model(
            model=test_model,
            graph="ASPASN-PWY",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
        )
        # Updating dummy biomass
        test_model.reactions.get_by_id("Biomass_c").add_metabolites(
            {
                test_model.metabolites.get_by_id("NICOTINE_c"): -0.25,
                test_model.metabolites.get_by_id("ASN_c"): -1,
            }
        )
        test_model.objective = "Biomass_c"
        self.assertGreater(abs(test_model.optimize().objective_value), 0)

    def test_cyclical_pathwways(self):
        test_model = main_model1.copy()
        # Adding Mannitol cycle
        ex.add_graph_to_model(
            model=test_model,
            graph="PWY-6531",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
        )
        # Adding glyoxylate cycle
        ex.add_graph_to_model(
            model=test_model,
            graph="GLYOXYLATE-BYPASS",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=[],
        )
        # Updating dummy biomass
        test_model.reactions.get_by_id("Biomass_c").add_metabolites(
            {
                test_model.metabolites.get_by_id("BETA_D_FRUCTOSE_c"): -1,
                test_model.metabolites.get_by_id("GLYOX_c"): -0.5,
            }
        )
        test_model.objective = "Biomass_c"
        self.assertGreater(abs(test_model.optimize().objective_value), 0)

    def test_multi_compartment(self):
        test_model = main_model1.copy()
        # This files has transports between cytosol and plastids
        add_reactions_from_file(
            model=test_model,
            filename=dir_input.joinpath("test_multi_reactions.txt"),
            database="META",
            directory=dir_data,
        )
        # ADDING COA as Exchange to omit its biosynthesis
        test_model.add_boundary(
            metabolite=test_model.metabolites.get_by_id("CO_A_e"),
            type="exchange",
        )
        ex.add_graph_to_model(
            model=test_model,
            graph="GLUTATHIONESYN-PWY",
            directory=dir_data,
            database="META",
            ignore_list=[],
            compartment="p",
        )
        test_model.reactions.get_by_id("Biomass_c").add_metabolites(
            {test_model.metabolites.get_by_id("GLUTATHIONE_p"): -1}
        )
        self.assertGreater(abs(test_model.optimize().objective_value), 0)
        _add_reaction_line_to_model(
            model=test_model,
            line="TRANS_GLY_cp, Transport GLY_cp | GLY_c:-1, GLY_p:1",
            directory=dir_data,
            database="META",
        )
        ex.add_graph_to_model(
            model=test_model,
            graph="PWY-1187",
            directory=dir_data,
            database="META",
            ignore_list=["PYRUVATE_c", "CO_A_c", "PROTON_c", "CPD_3746_c"],
            compartment="c",
        )
        test_model.reactions.get_by_id("Biomass_c").add_metabolites(
            {
                test_model.metabolites.get_by_id(
                    "3_METHYLSULFINYLPROPYL_GLUCOSINOLATE_c"
                ): -0.75,
                test_model.metabolites.get_by_id("GLUTATHIONE_p"): 1,
            }
        )
        self.assertGreater(abs(test_model.optimize().objective_value), 0)
        # Lipid initalization
        ex.add_graph_to_model(
            model=test_model,
            graph="PWY-4381",
            directory=dir_data,
            database="META",
            ignore_list=["CO_A_p"],
            compartment="p",
        )
        test_model.reactions.get_by_id("Biomass_c").add_metabolites(
            {test_model.metabolites.get_by_id("Acetoacetyl_ACPs_p"): -1}
        )
        # FIXME: get rid of extra sinks
        test_model.remove_reactions(["SK_Acetoacetyl_ACPs_p"])
        self.assertGreater(abs(test_model.optimize().objective_value), 0)


class TestingLargeModel(unittest.TestCase):
    def test_lineal_pathways(self):
        """In this test, unrelated pathways are appended to the large model
        and tested
        """
        test_model = main_model2.copy()
        # Autotroph enviroment
        test_model.reactions.Sucrose_tx.bounds = (0, 0)  # SUCROSE
        test_model.reactions.GLC_tx.bounds = (0, 0)  # GLUCOSE
        # # Photon uptake
        test_model.reactions.Photon_tx.bounds = (-1000, 250)
        # No ammonium
        test_model.reactions.NH4_tx.bounds = (0, 0)
        test_model.objective = "Biomass_tx"
        # In order to test, adding precursor as sink (glutathione)
        test_model.add_boundary(
            metabolite=test_model.metabolites.get_by_id("GLUTATHIONE_p"),
            type="sink",
        )
        # Adding methiin metabolism
        ex.add_graph_to_model(
            model=test_model,
            graph="PWY-7614",
            directory=dir_data,
            database="META",
            compartment="p",
            ignore_list=["GLUTATHIONE_p"],
        )
        # Checking demand for Methiin
        test_model.add_boundary(
            metabolite=test_model.metabolites.get_by_id("CPD_9277_p"),
            type="demand",
        )
        test_model.reactions.get_by_id("DM_CPD_9277_p").bounds = (1, 1000)
        self.assertGreater(test_model.optimize().fluxes["DM_CPD_9277_p"], 0)
        self.assertGreater(test_model.optimize().fluxes["RXN_8908_p"], 0)
        # Adding stachyose biosynthesis
        ex.add_graph_to_model(
            model=test_model,
            graph="PWY-5337",
            directory=dir_data,
            database="META",
            compartment="c",
            ignore_list=["PROTON_c", "MYO_INOSITOL_c"],
        )
        # Checking demand for Raffinose (created in the middle of pathway)
        test_model.add_boundary(
            metabolite=test_model.metabolites.get_by_id("CPD_1099_c"),
            type="demand",
        )
        test_model.reactions.get_by_id("DM_CPD_1099_c").bounds = (1, 1000)
        self.assertGreater(test_model.optimize().fluxes["DM_CPD_1099_c"], 0)
        self.assertGreater(test_model.optimize().fluxes["2.4.1.82_RXN_c"], 0)
        # Adding Abscisic acid
        test_model.add_boundary(
            metabolite=test_model.metabolites.get_by_id("CPD1F_133_p"),
            type="sink",
        )
        ex.add_graph_to_model(
            model=test_model,
            graph="PWY-695",
            directory=dir_data,
            database="META",
            compartment="p",
            ignore_list=["PROTON_p", "CPD1F_133_p"],
        )
        # FIXME: remove unneeded sinks
        test_model.remove_reactions(["SK_CPD_693_p"])
        # Checking demand for Methiin
        test_model.add_boundary(
            metabolite=test_model.metabolites.get_by_id("CPD_693_p"),
            type="demand",
        )
        test_model.reactions.get_by_id("DM_CPD_693_p").bounds = (1, 1000)
        self.assertGreater(test_model.optimize().fluxes["DM_CPD_693_p"], 0)
        # Reaction from middle section
        self.assertGreater(test_model.optimize().fluxes["RXN1F_155_p"], 0)
        # Last reaction
        self.assertGreater(test_model.optimize().fluxes["1.2.3.14_RXN_p"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
