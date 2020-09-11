#!/usr/bin/env python3
import cobra as cb
from pathlib import Path
from creation import add_meta_from_file, add_reactions_from_file, \
    _add_reaction_line_to_model
from pathways import add_graph_to_model

extraDir = Path.cwd().joinpath("extra")
test_dir = Path.cwd().joinpath("tests")
dirBiocyc = test_dir.joinpath("data").joinpath("biocyc")

# INFO:
# Model should have enough reactions/ metabolites for 2 short pathways
# No sink reactions, but instead a fake Biomass which get end metabolites.
# side-metabolites should be recycable
# Some metabolites such as ADP and ATP are being recycle by using redox
# reactions.

# New Model
test_model = cb.Model(
    id_or_model="test_model_v_1_0", name="Test model for cobraMod")
test_model.compartments = {
    "e": "extracellular",
    "c": "cytosol",
    "p": "plastid"}
add_meta_from_file(
    model=test_model, filename=extraDir.joinpath("metaToAdd_test_model.txt"),
    directory=dirBiocyc, database="META",)
add_reactions_from_file(
    model=test_model, filename=extraDir.joinpath(
        "reactionsToAdd_test_model.txt"),
    directory=dirBiocyc, database="META")
# Creating inital Exchanges reactions
for meta in test_model.metabolites:
    if meta.id[-1] == "e":
        test_model.add_boundary(metabolite=meta, type="exchange")
# Dummy objective reaction
test_model.objective = "EX_WATER_e"
# First pathway to create a proper dummy biomass reaction
# Nitrate reduction
add_graph_to_model(
    model=test_model, graph="PWY-381", directory=dirBiocyc,
    database="META",
    ignore_list=[])
# Creating new Dummy biomass reaction
_add_reaction_line_to_model(
    line="Biomass_c, Biomass reaction | GLN_c: -0.5",
    model=test_model, directory=dirBiocyc, database="META",)
test_model.objective = "Biomass_c"
test_model.objective_direction = "max"
# There is only one reaction. Thus, the sink needs to be removed
test_model.remove_reactions(["SK_GLN_c"])
# For first pathways in order to work
test_model.add_boundary(test_model.metabolites.get_by_id("CO_A_c"), "sink")
test_model.add_boundary(test_model.metabolites.get_by_id("CIT_c"), "sink")
for pathway in [
    "PWY-6964",  # Ammonium assimilation cycle
    "PWY-5690",  # TCA cycle
    "PYRUVDEHYD-PWY",
        ]:
    add_graph_to_model(
        model=test_model, graph=pathway, directory=dirBiocyc,
        database="META",
        ignore_list=["ETR_Quinols_c", "CIT_c"])
# Extending biomass
test_model.reactions.get_by_id("Biomass_c").add_metabolites({
            test_model.metabolites.get_by_id("GLT_c"): -1})
for pathway in ["CALVIN-PWY"]:
    add_graph_to_model(
        model=test_model, graph=pathway, directory=dirBiocyc,
        database="META",
        ignore_list=["GAP_c"])

for pathway in ["SUCSYN-PWY"]:
    add_graph_to_model(
        model=test_model, graph=pathway, directory=dirBiocyc,
        database="META",
        ignore_list=[
            "DIHYDROXY_ACETONE_PHOSPHATE_c",
            "PROTON_c",
            "GAPOXNPHOSPHN_RXN_c"])
# Extending biomass
test_model.reactions.get_by_id("Biomass_c").add_metabolites({
            test_model.metabolites.get_by_id("SUCROSE_c"): -1})
# For later tests
_add_reaction_line_to_model(
    line="UDPKIN-RXN, c",
    model=test_model, directory=dirBiocyc, database="META",)
test_model.remove_reactions(
    ["SK_UDP_c", "SK_UTP_c", "DM_SUCROSE_c"])
test_model.reactions.get_by_id("Biomass_c").bounds = (0.1, 1000)
# Saving model
cb.io.write_sbml_model(
    cobra_model=test_model,
    filename=str(
        test_dir.joinpath("input").joinpath("test_model01.sbml"))
)
