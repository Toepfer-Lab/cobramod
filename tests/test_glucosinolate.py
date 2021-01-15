from pathlib import Path
import unittest

import cobra as cb

from cobramod import extension as ex
from cobramod import creation as cr

dir_input = Path.cwd().joinpath("tests").joinpath("input")
path_model = dir_input.joinpath("test_model02.sbml")
dir_data = Path.cwd().joinpath("tests").joinpath("data")
dir_test_case = dir_input.joinpath("case_glucosinolate")
# Core Model
main_model = cb.io.read_sbml_model(filename=str(path_model))
test_model = main_model.copy()  # copy

if not dir_data.exists():
    dir_data.mkdir(parents=True)

# TODO: add option to stop at unbalanced reactions


def estimate_maintenance(light_intensity: float) -> float:
    """
    Random equation for ATPase constraint.
    """
    ATPase = (0.005 * light_intensity) + 2.8
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
            "R15-RXN": "RXN-15650",
        }
        self.ignore_list = [
            # no DM rxn or SK rxn will be created from these metabolites
            "WATER_c",
            "PROTON_c",
            "OXYGEN_MOLECULE_c",
            "CARBON_DIOXIDE_c",
            "3_5_ADP_c",
            "CPD_12607_c",
            "CPD_479_p",
            "2_Oxo_carboxylates_c",
            "Oxidized_hemoproteins_c",
        ]
        self.avoid_list = [
            # Reactions: This helps putting reactions into their correct
            # compartment.
            # GLS from tryptophan
            # oxidation in plastid
            "RXN_12061_c",
            "RXN_12062_c",
            "RXN_12063_c",
            # oxidation in ER
            "RXN_11413_c",  # From tryptophan
            "RXNQT_4311_c",  # Tetra
            "RXN_11417_c",  # Tetra
            "RXNQT_4312_c",  # Penta
            "RXN_11418_c",  # Penta
            "RXNQT_4313_c",  # Hexa
            "RXN_11419_c",
        ]  # Hexa

    def test_A_model_setting(self):
        # Autotroph enviroment
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
            directory=dir_data,
            database="ARA",
        )
        self.assertEqual(len(test_model.metabolites), 927)
        # Checking for Reactions
        self.assertEqual(len(test_model.reactions), 892)
        # This includes pathway from homophenylalanine
        cr.add_reactions_from_file(
            model=test_model,
            filename=dir_test_case.joinpath("new_reactions.txt"),
            database="ARA",
            directory=dir_data,
            replacement=self.replacement,
        )
        self.assertEqual(len(test_model.reactions), 959)
        # 3:1 rubisco rate
        Rubisco_rate = test_model.problem.Constraint(
            # Oxygenase
            3 * test_model.reactions.get_by_id("RXN_961_p").flux_expression
            # carboxylase
            - 1
            * test_model.reactions.get_by_id(
                "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"
            ).flux_expression,
            lb=0,
            ub=0,
        )  # to force it
        # Non-growth associated maintenance variables
        NGAM = test_model.problem.Constraint(
            # The sum of all 3 NADPH should be 1/3 of ATP
            3
            * (
                test_model.reactions.NADPHoxc_tx.flux_expression
                + test_model.reactions.NADPHoxm_tx.flux_expression
                + test_model.reactions.NADPHoxp_tx.flux_expression
            )
            - test_model.reactions.ATPase_tx.flux_expression,
            lb=0,
            ub=0,
        )
        # Defining ATPase constraint
        test_model.reactions.ATPase_tx.bounds = (
            estimate_maintenance(test_model.reactions.Photon_tx.upper_bound),
            estimate_maintenance(test_model.reactions.Photon_tx.upper_bound),
        )
        # Adding RuBisCO and NGAM constraint
        test_model.add_cons_vars([Rubisco_rate, NGAM])

    def test_B_precursors(self):
        # Glutathione synthesis
        ex.add_pathway(
            model=test_model,
            pathway="GLUTATHIONESYN-PWY",
            database="ARA",
            directory=dir_data,
            compartment="p",
            avoid_list=self.avoid_list,
        )
        # Original has an extra demand (total 962)
        self.assertEqual(len(test_model.reactions), 961)
        ex.add_pathway(
            model=test_model,
            pathway="PWY-5340",
            database="ARA",
            directory=dir_data,
            compartment="p",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        # One reaction was already in model
        self.assertEqual(len(test_model.reactions), 962)
        # Homomethionine synthesis
        ex.add_pathway(
            model=test_model,
            pathway="PWY-1186",
            database="ARA",
            directory=dir_data,
            compartment="p",
            replacement=self.replacement,
            avoid_list=["R15-RXN"],
            ignore_list=self.ignore_list,
        )
        # Its counterpart in cytosol was added from the file
        self.assertNotIn("R15_RXN_p", [rxn.id for rxn in test_model.reactions])
        # Methionine Elogation chain
        ex.add_pathway(
            model=test_model,
            pathway="PWYQT-4450",
            database="ARA",
            directory=dir_data,
            compartment="p",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertEqual(len(test_model.sinks), 0)

    def test_C_aliphatic(self):
        # From Homomethionine
        ex.add_pathway(
            model=test_model,
            pathway=[
                "RXN-11422",
                "RXN-11430",
                "RXN-11438",
                "RXN-2208",
                "RXN-2209",
                "RXN-2221",
            ],
            database="ARA",
            directory=dir_data,
            compartment="c",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertGreater(test_model.slim_optimize(), 0)
        # From Dihomemethionine
        ex.add_pathway(
            model=test_model,
            pathway=[
                "RXN-11423",
                "RXN-11431",
                "RXN-11439",
                "RXNQT-4324",
                "RXNQT-4329",
                "RXNQT-4334",
            ],
            database="ARA",
            directory=dir_data,
            compartment="c",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertGreater(test_model.slim_optimize(), 0)
        # From Trihomomethionine
        ex.add_pathway(
            model=test_model,
            pathway=[
                "RXN-11424",
                "RXN-11432",
                "RXN-11440",
                "RXNQT-4325",
                "RXNQT-4330",
                "RXNQT-4335",
            ],
            database="ARA",
            directory=dir_data,
            compartment="c",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertGreater(test_model.slim_optimize(), 0)
        # From Tetrahomomethionine
        ex.add_pathway(
            model=test_model,
            pathway="PWYQT-4473",
            database="ARA",
            directory=dir_data,
            compartment="c",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertGreater(test_model.slim_optimize(), 0)
        # From Pentahomomethionine
        ex.add_pathway(
            model=test_model,
            pathway="PWYQT-4474",
            database="ARA",
            directory=dir_data,
            compartment="c",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertGreater(test_model.slim_optimize(), 0)
        # From Hexahomomethionine
        ex.add_pathway(
            model=test_model,
            pathway="PWYQT-4475",
            database="ARA",
            directory=dir_data,
            compartment="c",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertGreater(test_model.slim_optimize(), 0)

    def test_D_indole(self):
        # Update metabolite
        meta = "S_ADENOSYLMETHIONINE_c"
        test_model.metabolites.get_by_id(meta).formula = "C15H23N6O5S"
        test_model.metabolites.get_by_id(meta).charge = 1
        # From Tryptophan
        ex.add_pathway(
            model=test_model,
            pathway="PWY-601",
            database="ARA",
            directory=dir_data,
            compartment="c",
            replacement=self.replacement,
            ignore_list=self.ignore_list,
            avoid_list=self.avoid_list,
        )
        self.assertGreater(test_model.slim_optimize(), 0)

    def test_E_biomass(self):
        test_biomass = cb.Reaction(id="Bio_opt_b", name="Bio_opt_b")
        # Biomass + Glucosinolates
        # BIOMASS = 3034.98641 + 17.92
        BIOMASS = 1
        with open(
            file=dir_test_case.joinpath("Bio_opt_rxn.txt"), mode="r"
        ) as f:
            for line in f:
                line = [part.strip().rstrip() for part in line.split(",")]
                if float(line[1]) > 0 or line[0] == "MALONYL_ACP_p":
                    continue
                else:
                    test_biomass.add_metabolites(
                        {
                            test_model.metabolites.get_by_id(line[0]): float(
                                line[1]
                            )
                            / BIOMASS
                        }
                    )
        test_model.add_reactions([test_biomass])
        self.assertGreater(test_model.slim_optimize(), 0)
        test_model.objective = "Bio_opt_b"
        test_fluxes = cb.flux_analysis.pfba(
            model=test_model, objective="Bio_opt_b"
        )
        self.assertAlmostEqual(
            test_fluxes.fluxes["Bio_opt_b"], 0.000897, delta=0.01
        )
        self.assertAlmostEqual(
            test_fluxes.fluxes["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"],
            16.981438106067703,
            delta=0.01,
        )
        self.assertAlmostEqual(
            test_fluxes.fluxes["RXN_961_p"], 5.660479368689235, delta=0.01
        )
        self.assertAlmostEqual(
            test_fluxes.fluxes["H2O_tx"], 13.016328186042248, delta=0.01
        )
        self.assertAlmostEqual(
            test_fluxes.fluxes["O2_tx"], -20.73902048985689, delta=0.01
        )
        self.assertAlmostEqual(
            test_fluxes.fluxes["SO4_tx"], 0.058057461417529495, delta=0.01
        )

    def test_E_glucosinolate(self):
        pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
