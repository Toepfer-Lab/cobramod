#!/usr/bin/env python3
from pathlib import Path
import unittest
# import pathways
import pathways as pt
# import xml.etree.ElementTree as ET
import cobra as cb
# from creation import get_xml_from_biocyc, add_reaction_from_root,\
#     add_meta_line_to_model, add_meta_from_file, add_reaction_line_to_model
# from itertools import chain

dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_biocyc = Path.cwd().joinpath("tests").joinpath("data").joinpath("biocyc")
main_model = cb.io.read_sbml_model(
    str(dir_input.joinpath("test_model01.sbml")))

if not dir_biocyc.exists():
    dir_biocyc.mkdir(parents=True)


class TestingModels(unittest.TestCase):

    def test_appeding_lineal_pathways(self):
        test_model = main_model.copy()
        test_model.add_boundary(test_model.metabolites.get_by_id(
            "ATP_c"), "sink")
        pt.add_graph_from_root(
            model=test_model, root="GLUCONEO-PWY", directory=dir_biocyc,
            ignore_list=["PYRUVATE_c", "MAL_c", "PROTON_c", "Pi_c", "ATP_c"])
        test_model.remove_reactions(["SK_ATP_c"])
        pass


if __name__ == "__main__":
    unittest.main(verbosity=2)