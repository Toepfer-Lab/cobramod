from cobramod import add_reactions_from_file
from cobra import Model
from cobra.test import create_test_model
from pathlib import Path

dir_data = Path.cwd().joinpath("data")
filename = Path.cwd().joinpath("bigg_to_meta.txt")

test_model = Model(id_or_model="test")
add_reactions_from_file(
    model=test_model, filename=filename, directory=dir_data, database="META"
)

test = create_test_model("textbook")
for rxn in test.reactions:
    print(rxn.id)

