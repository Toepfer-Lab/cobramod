from pathlib import Path

from cobra.io import read_sbml_model
from cobra.test import create_test_model

data_dir = Path(__file__).parents[0].joinpath("data")

textbook = create_test_model(model_name="textbook")

textbook_biocyc = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_biocyc.sbml"))
)
textbook_kegg = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_kegg.sbml"))
)
