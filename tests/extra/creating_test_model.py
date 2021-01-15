#!/usr/bin/env python3
import cobra as cb
from pathlib import Path
from cobramod.extension import add_pathway
from cobramod.creation import (
    add_meta_from_file,
    add_reactions_from_file,
    _add_reaction_line_to_model,
)

test_dir = Path.cwd().joinpath("tests")
extra_dir = test_dir.joinpath("extra")
dir_data = test_dir.joinpath("data")

if not dir_data.exists():
    dir_data.mkdir(parents=True)

# INFO:
# Model should have enough reactions/ metabolites for 2 short pathways
# No sink reactions, but instead a fake Biomass which get end metabolites.
# side-metabolites should be recycable
# Some metabolites such as ADP and ATP are being recycle by using redox
# reactions.

# New Model
test_model = cb.Model(
    id_or_model="test_model_v_1_0", name="Test model for cobraMod"
)
test_model.compartments = {
    "e": "extracellular",
    "c": "cytosol",
    "p": "plastid",
}
add_meta_from_file(
    model=test_model,
    filename=extra_dir.joinpath("metaToAdd_test_model.txt"),
    directory=dir_data,
    database="META",
)
add_reactions_from_file(
    model=test_model,
    filename=extra_dir.joinpath("reactionsToAdd_test_model.txt"),
    directory=dir_data,
    database="META",
)
# Creating inital Exchanges reactions
for meta in test_model.metabolites:
    if meta.id[-1] == "e":
        test_model.add_boundary(metabolite=meta, type="exchange")
# Dummy objective reaction
test_model.objective = "EX_WATER_e"
# First pathway to create a proper dummy biomass reaction
# Nitrate reduction
add_pathway(
    model=test_model,
    pathway="PWY-381",
    directory=dir_data,
    database="META",
    compartment="c",
    ignore_list=[],
)
# Creating new Dummy biomass reaction
_add_reaction_line_to_model(
    line="Biomass_c, Biomass reaction | GLN_c: -0.5",
    model=test_model,
    directory=dir_data,
    database="META",
)
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
    add_pathway(
        model=test_model,
        pathway=pathway,
        directory=dir_data,
        database="META",
        compartment="c",
        ignore_list=["ETR_Quinols_c", "CIT_c"],
    )
# Extending biomass
test_model.reactions.get_by_id("Biomass_c").add_metabolites(
    {test_model.metabolites.get_by_id("GLT_c"): -1}
)
for pathway in ["CALVIN-PWY"]:
    add_pathway(
        model=test_model,
        pathway=pathway,
        directory=dir_data,
        database="META",
        compartment="c",
        ignore_list=["GAP_c"],
    )

for pathway in ["SUCSYN-PWY"]:
    add_pathway(
        model=test_model,
        pathway=pathway,
        directory=dir_data,
        database="META",
        compartment="c",
        ignore_list=[
            "DIHYDROXY_ACETONE_PHOSPHATE_c",
            "PROTON_c",
            "GAPOXNPHOSPHN_RXN_c",
        ],
    )
# Extending biomass
test_model.reactions.get_by_id("Biomass_c").add_metabolites(
    {test_model.metabolites.get_by_id("SUCROSE_c"): -1}
)
# For later tests
_add_reaction_line_to_model(
    line="UDPKIN-RXN, c", model=test_model, directory=dir_data, database="META"
)
test_model.remove_reactions(["SK_UDP_c", "SK_UTP_c"])
test_model.reactions.get_by_id("Biomass_c").bounds = (0.1, 1000)
# Saving model
cb.io.write_sbml_model(
    cobra_model=test_model,
    filename=str(test_dir.joinpath("input").joinpath("test_model01.sbml")),
)
