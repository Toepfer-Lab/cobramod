import pathways as pt
import creation as cr
import unittest
from pathlib import Path
import cobra as cb

dir_input = Path.cwd().joinpath("tests").joinpath("input")
path_model = dir_input.joinpath("test_model02.sbml")
dir_biocyc = Path.cwd().joinpath("tests").joinpath("data").joinpath("biocyc")
dir_test_case = dir_input.joinpath("case_glucosinolate")
# Core Model
main_model = cb.io.read_sbml_model(filename=str(path_model))


def estimate_maintenance(light_intensity: float) -> float:
    """Equation for ATPase constraint.

    :param light_intensity: maximum uptake of photons (in mikromol/m^2s)
    :type light_intensity: float or int
    :return: Energy constraint for ATPase.
    :rtype: float
    """
    ATPase = (0.0049 * light_intensity) + 2.7851
    return float(ATPase)


class GlucosinolateTest(unittest.TestCase):

    def setUp(self):
        self.replacement = {
            # Metabolites
            "CPD-12575": "UDP-GLUCOSE",
            "Red-NADPH-Hemoprotein-Reductases": "Reduced-hemoproteins",
            "Ox-NADPH-Hemoprotein-Reductases": "Oxidized-hemoproteins",
            # All generic reactions showed Glutamine as commom aminoacid
            "Amino-Acids-20": "GLT",
            "2-Oxo-carboxylates": "2-KETOGLUTARATE",
            # Reactions (generic to specific)
            "R15-RXN": "RXN-15650"}
        self.ignore_list = [
            # no DM rxn or SK rxn will be created from these metabolites
            "WATER_c", "PROTON_c", "OXYGEN_MOLECULE_c",
            "CARBON_DIOXIDE_c", "3_5_ADP_c", "CPD_12607_c", "CPD_479_p",
            "2_Oxo_carboxylates_c"]
        self.avoid_list = [
            # Reactions: This helps putting reactions into their correct
            # compartment.
            # GLS from tryptophan
            # oxidation in plastid
            "RXN_12061_c", "RXN_12062_c", "RXN_12063_c",
            "RXN_11413_c"  # oxidation in ER
            ]
        self.test_model = main_model.copy()  # copy
        test_model = self.test_model
        # objective = "Biomass_tx"
        # SETTING ENVIROMENT
        # Autotroph
        test_model.reactions.Sucrose_tx.bounds = (0, 0)  # SUCROSE
        test_model.reactions.GLC_tx.bounds = (0, 0)  # GLUCOSE
        # # Photon uptake
        test_model.reactions.Photon_tx.bounds = (-1000, 250)
        # No ammonium
        test_model.reactions.NH4_tx.bounds = (0, 0)
        # Checking for Metabolites
        self.assertEqual(len(test_model.metabolites), 861)
        cr.add_meta_from_file(
            model=test_model,
            filename=dir_test_case.joinpath("new_metabolites.txt"),
            directory=dir_biocyc,
            database="ARA")
        self.assertEqual(len(test_model.metabolites), 927)
        # Checking for Reactions
        self.assertEqual(len(test_model.reactions), 892)
        cr.add_reaction_from_file(
            model=test_model,
            filename=dir_test_case.joinpath("new_reactions.txt"),
            database="ARA",
            directory=dir_biocyc,
            replacement_dict=self.replacement)
        self.assertEqual(len(test_model.reactions), 959)
        # 3:1 rubisco rate
        Rubisco_rate = test_model.problem.Constraint(
            3 * test_model.reactions.RXN_961_p.flux_expression -   # Oxygenase
            1 * test_model.reactions.\
            # carboxylase
            RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p.flux_expression,
            lb=0, ub=0)  # to force it
        # Non-growth associated maintenance variables
        NGAM = test_model.problem.Constraint(
            # The sum of all 3 NADPH should be 1/3 of ATP
            3*(test_model.reactions.NADPHoxc_tx.flux_expression +
                test_model.reactions.NADPHoxm_tx.flux_expression +
                test_model.reactions.NADPHoxp_tx.flux_expression) -
            test_model.reactions.ATPase_tx.flux_expression,
            lb=0, ub=0)
        # Defining ATPase constraint
        test_model.reactions.ATPase_tx.bounds = (
            estimate_maintenance(test_model.reactions.Photon_tx.upper_bound),
            estimate_maintenance(test_model.reactions.Photon_tx.upper_bound))
        # Adding RuBisCO and NGAM constraint
        test_model.add_cons_vars([
            Rubisco_rate, NGAM])
        # self.test_model = test_model

    def test_A_precursors(self):
        replacement = self.replacement
        ignore_list = self.ignore_list
        # Adding new pathways
        # test_model = self.test_model
        # Glutathione synthesis
        pt.add_graph_to_model(
            model=self.test_model, graph="GLUTATHIONESYN-PWY", database="ARA",
            directory=dir_biocyc, compartment="p")
        # Original has an extra demand (total 962)
        self.assertEqual(len(self.test_model.reactions), 961)
        pt.add_graph_to_model(
            model=self.test_model, graph="PWY-5340", database="ARA",
            directory=dir_biocyc, compartment="p",
            replacement_dict=replacement, ignore_list=ignore_list)
        # One reaction was already in model
        self.assertEqual(len(self.test_model.reactions), 962)
        # Homomethionine synthesis
        pt.add_graph_to_model(
            model=self.test_model, graph="PWY-1186", database="ARA",
            directory=dir_biocyc, compartment="p",
            replacement_dict=replacement,
            avoid_list=["R15-RXN"], ignore_list=ignore_list)
        # Its counterpart in cytosol was added from the file
        self.assertNotIn(
            "R15_RXN_p",
            [rxn.id for rxn in self.test_model.reactions])
        # Methionine Elogation chain
        pt.add_graph_to_model(
            model=self.test_model, graph="PWYQT-4450", database="ARA",
            directory=dir_biocyc, compartment="p",
            replacement_dict=replacement, ignore_list=ignore_list)
        self.assertEqual(len(self.test_model.sinks), 0)
        # self.test_model = test_model

    def test_B_aliphatic(self):
        replacement = self.replacement
        ignore_list = self.ignore_list
        # Adding new pathways
        # test_model = self.test_model
        # From Homomethionine
        # pt.add_graph_to_model(
        #     model=test_model, graph=[
        #         'RXN-11422', 'RXN-11430', 'RXN-11438', 'RXN-2208',
        #         'RXN-2209', 'RXN-2221'],
        #     database="ARA",
        #     directory=dir_biocyc, compartment="c",
        #     replacement_dict=replacement, ignore_list=ignore_list)
        # # From Dihomemethionine
        # pt.add_graph_to_model(
        #     model=test_model, graph=[
        #         'RXN-11423', 'RXN-11431',
        #         'RXN-11439', 'RXNQT-4324', 'RXNQT-4329', 'RXNQT-4334'],
        #     database="ARA",
        #     directory=dir_biocyc, compartment="c",
        #     replacement_dict=replacement, ignore_list=ignore_list)
        # # From Trihomomethionine
        # pt.add_graph_to_model(
        #     model=test_model, graph=[
        #         'RXN-11424', 'RXN-11432', 'RXN-11440',
        #         'RXNQT-4325', 'RXNQT-4330', 'RXNQT-4335'],
        #     database="ARA",
        #     directory=dir_biocyc, compartment="c",
        #     replacement_dict=replacement, ignore_list=ignore_list)
        # # From Tetrahomomethionine
        # pt.add_graph_to_model(
        #     model=test_model, graph="PWYQT-4473",
        #     database="ARA",
        #     directory=dir_biocyc, compartment="c",
        #     replacement_dict=replacement, ignore_list=ignore_list)
        # # From Pentahomomethionine
        # pt.add_graph_to_model(
        #     model=test_model, graph="PWYQT-4474",
        #     database="ARA",
        #     directory=dir_biocyc, compartment="c",
        #     replacement_dict=replacement, ignore_list=ignore_list)
        # # From Hexahomomethionine
        # pt.add_graph_to_model(
        #     model=test_model, graph="PWYQT-4475",
        #     database="ARA",
        #     directory=dir_biocyc, compartment="c",
        #     replacement_dict=replacement, ignore_list=ignore_list)

        # pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
