#!/usr/bin/env python3
from logging import DEBUG
from pathlib import Path
import unittest

from cobra import Metabolite, Reaction

import cobramod.mod_parser as md
from cobramod.debug import debug_log

debug_log.setLevel(DEBUG)
dir_data = Path.cwd().joinpath("tests").joinpath("data")


class RetrievalTesting(unittest.TestCase):
    def test_get_data(self):
        # CASE 1a: Simple retrieval from Biocyc (metabolite)
        test_dict = md.get_data(
            directory=dir_data, database="META", identifier="WATER"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 1b: simple retrieval from Biocyc (reaction)
        test_dict = md.get_data(
            directory=dir_data, database="META", identifier="GLYOXII-RXN"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        # CASE 2: simple retrieval from AraCyc
        # TODO: Case 2
        # CASE 3a: Simple retrieval from KEGG (metabolite)
        test_dict = md.get_data(
            directory=dir_data, database="KEGG", identifier="C00001"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 3b: Simple retrieval from KEGG (Reaction)
        test_dict = md.get_data(
            directory=dir_data, database="KEGG", identifier="R02736"
        )
        # CASE 3c: Simple retrieval from KEGG (Pathway)
        test_dict = md.get_data(
            directory=dir_data, database="KEGG", identifier="M00001"
        )
        self.assertEqual(first=test_dict["TYPE"], second="Pathway")
        self.assertEqual(first=test_dict["ENTRY"], second="M00001")
        # CASE 4a: Simple retrieval from BIGG (metabolite)
        test_dict = md.get_data(
            directory=dir_data,
            database="BIGG",
            identifier="accoa_c",
            model_id="e_coli_core",
        )
        # CASE 4b: Simple retrieval from BIGG (reaction)
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        self.assertEqual(first=test_dict["ENTRY"], second="accoa")
        test_dict = md.get_data(
            directory=dir_data,
            database="BIGG",
            identifier="CS",
            model_id="e_coli_core",
        )
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        self.assertEqual(first=test_dict["ENTRY"], second="CS")

    def test_check_attribute(self):
        # CASE: Generic compounds:
        for test_dict in (
            md.get_data(
                directory=dir_data,
                database="META",
                identifier="Polyphosphates",
            ),
            md.get_data(
                directory=dir_data, database="KEGG", identifier="C00404"
            ),
            md.get_data(
                directory=dir_data, database="KEGG", identifier="C00139"
            ),
        ):
            self.assertEqual(first=test_dict["FORMULA"], second="X")

    def test_translate(self):
        # CASE 1: Regular compound KEGG
        test_string = md.translate(
            directory=dir_data, target="C00002", database="CAS"
        )
        self.assertEqual(first=test_string, second="56-65-5")
        # CASE 1: Regular compound Biocyc
        test_string = md.translate(
            directory=dir_data, target="AMP", database="PUBCHEM"
        )
        self.assertEqual(first=test_string, second="15938965")
        # CASE 1: Regular compound BIGG
        test_string = md.translate(
            directory=dir_data, target="accoa_c", database="CHEBI"
        )
        self.assertEqual(first=test_string, second="CHEBI:13712")

    def test__valid_translation(self):
        # CASE 1a: Checking water from Biocyc to KEGG
        test_metabolite = Metabolite(id="C00001_c")
        self.assertTrue(
            expr=md._valid_translation(
                metabolite=test_metabolite, pattern="WATER", directory=dir_data
            )
        )
        # CASE 1b: Checking water from BIGG to KEGG
        test_metabolite = Metabolite(id="C00001_c")
        self.assertTrue(
            expr=md._valid_translation(
                metabolite=test_metabolite, pattern="h2o", directory=dir_data
            )
        )
        # CASE 1c: Checking water from biocyc to unknown database
        test_metabolite = Metabolite(id="7732_18_5_c")
        self.assertTrue(
            expr=md._valid_translation(
                metabolite=test_metabolite, pattern="WATER", directory=dir_data
            )
        )
        # CASE 2a: Checking Reaction from BIGG to KEGG
        test_metabolite = Reaction(id="R00228_c")
        self.assertTrue(
            expr=md._valid_translation(
                metabolite=test_metabolite, pattern="ACALD", directory=dir_data
            )
        )
        # CASE 2b: Checking Reaction from Biocyc to KEGG
        test_metabolite = Reaction(id="R00228_c")
        self.assertTrue(
            expr=md._valid_translation(
                metabolite=test_metabolite,
                pattern="ACETALD-DEHYDROG-RXN",
                directory=dir_data,
            )
        )
        # CASE 2c: Checking Reaction from  Unknown (RHEA) to BIGG
        test_metabolite = Reaction(id="23288_c")
        self.assertTrue(
            expr=md._valid_translation(
                metabolite=test_metabolite, pattern="ACALD", directory=dir_data
            )
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
