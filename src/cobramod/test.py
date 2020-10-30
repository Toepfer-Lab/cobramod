from pathlib import Path
from cobra.io import read_sbml_model

data_dir = Path.cwd().joinpath("src").joinpath("cobramod").joinpath("data")

mini_model = read_sbml_model(
    filename=str(data_dir.joinpath("mini_biocyc.sbml"))
)
