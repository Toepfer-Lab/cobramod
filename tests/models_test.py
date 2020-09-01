#!/usr/bin/env python3
from pathlib import Path
import unittest
# import pathways
import pathways as pt
# import xml.etree.ElementTree as ET
import cobra as cb
from creation import add_reaction_line_to_model

dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_biocyc = Path.cwd().joinpath("tests").joinpath("data").joinpath("biocyc")
main_model = cb.io.read_sbml_model(
    str(dir_input.joinpath("test_model01.sbml")))

if not dir_biocyc.exists():
    dir_biocyc.mkdir(parents=True)


class TestingModels(unittest.TestCase):

    def test_appending_lineal_pathways(self):
        """For this test, Gluconeogenesis; L-aspartate and L-asparagine
        biosynthesis, and Nicotine biosynthesis added to test for feasible
        results. Due to a small model,  ATP needs to enter as a sink reaction
        but is later removed from model without interfering with the the new
        pathways
        """
        test_model = main_model.copy()
        test_model.add_boundary(test_model.metabolites.get_by_id(
            "ATP_c"), "sink")
        # Adding Gluconeogenesis
        pt.add_graph_from_root(
            model=test_model, root="GLUCONEO-PWY", directory=dir_biocyc,
            ignore_list=["MAL_c", "PROTON_c", "Pi_c", "ATP_c"])
        test_model.remove_reactions([
            "SK_ATP_c"])
        self.assertGreater(
            abs(test_model.optimize().objective_value), 0)
        add_reaction_line_to_model(
            model=test_model,
            line=(
                "Redox_Hemoprotein_c, Redox_Hemoprotein_c | " +
                "Ox-NADPH-Hemoprotein-Reductases_c:-1, " +
                "Red-NADPH-Hemoprotein-Reductases_c: 1"),
            directory=dir_biocyc)
        # Adding Nicotine pathway
        pt.add_graph_from_root(
            model=test_model, root="PWY-5316", directory=dir_biocyc,
            ignore_list=[])
        # Ading aspartate and asparagine biosynthesis
        pt.add_graph_from_root(
            model=test_model, root="ASPASN-PWY", directory=dir_biocyc,
            ignore_list=[])
        # Updating dummy biomass
        test_model.reactions.get_by_id("Biomass_c").add_metabolites({
            test_model.metabolites.get_by_id("NICOTINE_c"): -0.25,
            test_model.metabolites.get_by_id("ASN_c"): -1})
        test_model.objective = "Biomass_c"
        self.assertGreater(
            abs(test_model.optimize().objective_value), 0)

    def test_cyclical_pathwways(self):
        test_model = main_model.copy()
        # Adding Mannitol cycle
        pt.add_graph_from_root(
            model=test_model, root="PWY-6531", directory=dir_biocyc,
            ignore_list=[])
        # Adding glyoxylate cycle
        pt.add_graph_from_root(
            model=test_model, root="GLYOXYLATE-BYPASS", directory=dir_biocyc,
            ignore_list=[])
        # Updating dummy biomass
        test_model.reactions.get_by_id("Biomass_c").add_metabolites({
            test_model.metabolites.get_by_id("MANNITOL_c"): -0.25,
            test_model.metabolites.get_by_id("GLYOX_c"): -0.5
            })
        test_model.objective = "Biomass_c"
        self.assertGreater(
            abs(test_model.optimize().objective_value), 0)

    def test_multi_compartment(self):
        test_model = main_model.copy()
        # add_meta_line_to_model(
        #     line="ACP, c", model=test_model, directory=dir_biocyc)
        test_model.add_boundary(test_model.metabolites.get_by_id(
            "CO_A_c"), "sink")
        add_reaction_line_to_model(
            model=test_model,
            line=(
                "Redox_Hemoprotein_c, Redox_Hemoprotein_c | " +
                "Ox-NADPH-Hemoprotein-Reductases_c:-1, " +
                "Red-NADPH-Hemoprotein-Reductases_c: 1"),
            directory=dir_biocyc)
        pt.add_graph_from_root(
            model=test_model, root="GLUTATHIONESYN-PWY", directory=dir_biocyc,
            ignore_list=["PYRUVATE_c"], compartment="c")
        # FIXME: remove demands
        test_model.remove_reactions(["DM_GLUTATHIONE_c"])
        pt.add_graph_from_root(
            model=test_model, root="PWY-1187", directory=dir_biocyc,
            ignore_list=["PYRUVATE_c", "CO_A_c"], compartment="c")
        pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
