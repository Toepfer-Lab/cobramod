import webbrowser
from pathlib import Path
from random import randint
from unittest import TestCase

from cobra import Reaction, Metabolite

from cobramod import add_pathway
from cobramod.test import textbook_biocyc
from cobramod.visualization.threejs import ForceGraphBuilder

dir_data = Path(__file__).resolve().parent.joinpath("data")
dir_input = Path(__file__).resolve().parent.joinpath("input")


class TestForceGraphBuilder(TestCase):

    def test_pathway(self):
        test_model = textbook_biocyc.copy()
        add_pathway(
            model=test_model,
            pathway="SALVADEHYPOX-PWY",
            compartment="c",
            directory=dir_data,
            database="ECO",
            ignore_list=[],
            show_imbalance=False,
        )
        # Test fluxes
        test_pathway = test_model.groups.get_by_id("SALVADEHYPOX-PWY")

        builder = ForceGraphBuilder(test_pathway)
        builder.open_html()

        test_solution = {
            reaction.id: randint(-4, 4) for reaction in test_pathway.members
        }

        builder = ForceGraphBuilder(test_pathway, solution= test_solution)
        builder.open_html()