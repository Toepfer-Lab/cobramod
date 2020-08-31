#!/usr/bin/env python3
import cobra as cb
from pathlib import Path
import logging
from creation import add_meta_from_file, add_reaction_from_file, \
    add_reaction_line_to_model
from pathways import add_graph_from_root

# Creating corresponding Logs
# Format
DebugFormatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# Handler
DebugHandler = logging.FileHandler("debug.log", mode="a+")  
DebugHandler.setFormatter(DebugFormatter)
# Log
DebugLog = logging.getLogger("DebugLog")
DebugLog.setLevel(logging.DEBUG)
# Directories
extraDir = Path.cwd().joinpath("extra")
dirBiocyc = Path.cwd().joinpath("tests").joinpath("data").joinpath("biocyc")

# INFO: !!
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
    directory=dirBiocyc)
add_reaction_from_file(
    model=test_model, filename=extraDir.joinpath(
        "reactionsToAdd_test_model.txt"),
    directory=dirBiocyc)
# Creating inital Exchanges reactions
for meta in test_model.metabolites:
    if meta.id[-1] == "e":
        test_model.add_boundary(metabolite=meta, type="exchange")
# For first pathways in order to work
test_model.add_boundary(test_model.metabolites.get_by_id("MAL_c"), "sink")
test_model.add_boundary(test_model.metabolites.get_by_id("CO_A_c"), "sink")
# test_model.add_boundary(test_model.metabolites.get_by_id("GAP_c"), "sink")
for pathway in [
    "PWY-381",  # Nitrate reduction
    "PWY-6964",  # Ammonium assimilation cycle
    "PWY-5690",  # TCA cycle
    # "CALVIN-PWY"  # Calvin cycle
    "PYRUVDEHYD-PWY",
    # "PWY-5172",
    ]:
# for pathway in []:
    add_graph_from_root(
        model=test_model, root=pathway, directory=dirBiocyc,
        ignore_list=[
            "MAL_c", "CO_A_c", "ETR_Quinols_c",
            "", "GAP_c"
            ])
add_reaction_line_to_model(
    line="Biomass_c, Biomass reaction | ATP_c:-0.5",
    model=test_model)
test_model.remove_reactions(
    ["SK_CO_A_c", "SK_MAL_c"])
test_model.reactions.get_by_id("Biomass_c").bounds = (0.1, 1000)
test_model.objective = "Biomass_c"
test_model.objective_direction = "max"
# # Saving model
cb.io.write_sbml_model(
    cobra_model=test_model,
    filename=str(extraDir.joinpath("test_model.sbml"))
)
