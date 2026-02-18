import logging
from unittest import TestCase

import cobra.io
from cobra import Metabolite, Model

from cobramod.dbwalker.crossreferences import add_crossreferences2metabolite, add_crossreferences, id2annotation

#logging.basicConfig(
#    level=logging.DEBUG,
#    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
#    handlers=[logging.StreamHandler()],
#)

logger = logging.getLogger("cobramod")
logger.setLevel(logging.DEBUG)


class TestCrossReferences(TestCase):

    def test_example_model_with5Metabolite(self):
        model = Model()

        metabolite = Metabolite(id = "BioCyc")
        metabolite.annotation = {"biocyc": "ECOLI:CPD-12457"}
        model.add_metabolites([metabolite])

        metabolite = Metabolite(id = "CHEBI")
        metabolite.annotation = {"CHEBI": "16454"}
        model.add_metabolites([metabolite])

        metabolite = Metabolite(id="PubChem")
        metabolite.annotation = {"pubchem.compound": "24578"}
        model.add_metabolites([metabolite])

        metabolite = Metabolite(id="Kegg")
        metabolite.annotation = {"kegg.compound": "C12458"}
        model.add_metabolites([metabolite])


        add_crossreferences(model, save_log= "./annotation_log.csv")

        print(metabolite.annotation)
        self.assertTrue(True)

    def test_add_crossreferences2metabolite(self):
        """Test adding cross-references to a metabolite."""

        cid = "2519"

        metabolite = Metabolite()
        metabolite.id = cid
        metabolite.annotation = {"biocyc": "ECO:CYTOSINE"}

        add_crossreferences2metabolite(metabolite)

        print(metabolite.annotation)
        self.assertTrue(True)

    def test_add_crossreferences(self):

        model = cobra.io.read_sbml_model("/home/jan/Work/IPK/cobramod2/tests/dbwalker/model_Pmegm12_B_E_Growth185.xml")


        #print(model.metabolites[0])
        #print(model.metabolites[0].annotation)

        #id2annotation(model)

        #print(model.metabolites[0].annotation)

        #add_crossreferences2metabolite(model.metabolites[0])

        #print(model.metabolites[0].annotation)

        add_crossreferences(model, save_log= "./annotation_log.csv")

        cobra.io.write_sbml_model("/home/jan/Work/IPK/cobramod2/tests/dbwalker/model_Pmegm12_B_E_Growth185_annotation.xml")

    def test_another(self):
        model = cobra.io.read_sbml_model("/home/jan/Work/IPK/cobramod2/tests/dbwalker/model_Pmegm12_B_E_Growth185.xml")

        print(len(model.metabolites))