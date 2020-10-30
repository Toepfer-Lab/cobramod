#!/usr/bin/env python3
# from cobramod.parsing.bigg import _get_json_bigg, _get_db_names
from cobra.test import create_test_model
from cobra.io import write_sbml_model
from cobramod.utils import create_replacement
from pathlib import Path

dir_data = Path.cwd().joinpath("data")
filename = Path.cwd().joinpath("bigg_to_meta.txt")

test = create_test_model("textbook")

# with open(file=str(filename), mode="w+") as f:
#     for meta in test.metabolites:
#         try:
#             data_dict = _get_json_bigg(
#                 directory=dir_data,
#                 identifier=meta.id,
#                 model_id="e_coli_core",
#                 object_type="metabolite",
#             )
#             replacement = _get_db_names(data_dict=data_dict, database="META")
#             replacement = (
#                 replacement[5:].replace("-", "_") + meta.id[-2:]
#             )
#         except KeyError:
#             replacement = "Not found"
#         msg = f"{meta.id}: {replacement}\n"
#         f.write(msg)


test_dict = create_replacement(filename=filename)
for meta in test.metabolites:
    meta.id = test_dict[meta.id]
    print(meta.id)

write_sbml_model(cobra_model=test, filename="mini_biocyc.sbml")
