from cobra import Model
from cobra.io import write_sbml_model
from cobramod import add_reactions_from_file
from pathlib import Path

dir_test = Path.cwd().joinpath("tests")
dir_data = dir_test.joinpath("data")
file_path = dir_test.joinpath("extra").joinpath("rxn_mini.txt")

test_model = Model(
    id_or_model="mini_biocyc", name="COBRApy Mini converted to Biocyc"
)
test_model.compartments = {"e": "extracellular", "c": "cytosol"}
add_reactions_from_file(
    model=test_model, filename=file_path, directory=dir_data, database="META"
)
write_sbml_model(
    cobra_model=test_model,
    filename=str(dir_test.joinpath("extra").joinpath("mini_biocyc.sbml")),
)
