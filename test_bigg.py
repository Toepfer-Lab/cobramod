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

for meta in test_model.reactions.Biomass_c.metabolites:
    if len(meta.reactions) == 1:
        print(meta.id)
        print(meta.reactions)
test_model.objective = "Biomass_c"
test_model.objective_direction = "max"
print(test_model.optimize())
print("test")
