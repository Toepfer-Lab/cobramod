import logging
from unittest import TestCase

import cobra.io
from cobra import Metabolite, Model

from cobramod.dbwalker.crossreferences import (
    add_crossreferences,
    add_crossreferences2metabolite,
    id2annotation,
)

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

        metabolite_kegg = Metabolite(id="Kegg")
        metabolite_kegg.annotation = {"kegg.compound": "C12458"}
        model.add_metabolites([metabolite_kegg])


        add_crossreferences(model, save_log= "./annotation_log.csv")

        print(metabolite_kegg.annotation)
        metabolite_kegg_expected_annotation = {
            "kegg.compound": "C12458",
            "pubchem.compound": "443852",
            "inchi": "InChI=1S/C12H14O3/c1-8(2)3-4-9-7-10(12(14)15)5-6-11(9)13/h3,5-7,13H,4H2,1-2H3,(H,14,15)",
            "inchikey": "LBSJJNAMGVDGCU-UHFFFAOYSA-N",
            "smiles": "CC(=CCC1=C(C=CC(=C1)C(=O)O)O)C",
        }

        self.assertDictEqual(metabolite_kegg.annotation, metabolite_kegg_expected_annotation)

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