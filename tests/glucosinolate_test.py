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
        test_model = main_model.copy()  # copy
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
            directory=dir_biocyc)
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
        self.test_model = test_model

    def test_yeh(self):
        # Adding new pathways
        replacement = {
            "R15-RXN": "RXN-15650"
        }
        test_model = self.test_model
        # Precursors
        pt.add_graph_from_root(
            model=test_model, root="GLUTATHIONESYN-PWY", database="ARA",
            directory=dir_biocyc, compartment="p")
        # Original has an extra demand (total 962)
        self.assertEqual(len(test_model.reactions), 961)
        pt.add_graph_from_root(
            model=test_model, root="PWY-5340", database="ARA",
            directory=dir_biocyc, compartment="p")
        # One reaction was already in model
        self.assertEqual(len(test_model.reactions), 962)
        pt.add_graph_from_root(
            model=test_model, root="PWY-1186", database="ARA",
            directory=dir_biocyc, compartment="p",
            replacement_dict=replacement)
        self.assertIn(
            "RXN_15650_p",
            [rxn.id for rxn in test_model.reactions])
        # Test
        # test_model.objective = objective
        # test_solution = cb.flux_analysis.pfba(
        #     test_model, objective=objective)
        # test_model.summary(
        #     solution=test_solution)
        pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
