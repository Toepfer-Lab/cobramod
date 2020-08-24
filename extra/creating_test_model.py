#!/usr/bin/env python3
import cobra as cb
from pathlib import Path
import logging
from GLS import add_meta_from_file, add_reaction_from_file
from pathways import testAndAddCompletePathway

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
dirBiocyc = Path.cwd().joinpath("biocyc")

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
test_model.add_boundary(test_model.metabolites.get_by_id("MAL_c"), "sink")
test_model.add_boundary(test_model.metabolites.get_by_id("CO_A_c"), "sink")
for pathway in [
    # "PWY-381",
    "CALVIN-PWY"
    # "PWY-6964",
    # "PWY-5690",
    # "PWY-5172",
    # "PYRUVDEHYD-PWY"
    ]:
# for pathway in []:
    testAndAddCompletePathway(
        model=test_model, xmlRoot=pathway, directory=dirBiocyc,
        ignoreList=[
            "MAL_c", "CO_A_c", "ETR_Quinols_c",
            "Pyruvate_dehydrogenase_dihydrolipoate_c"
            ])
test_model.reactions.get_by_id("Biomass_c").bounds = (0.1, 1000)
test_model.objective_direction = "max"
test_model.objective = "Biomass_c"
# # Saving model
cb.io.write_sbml_model(
    cobra_model=test_model,
    filename=str(extraDir.joinpath("test_model.sbml"))
)
