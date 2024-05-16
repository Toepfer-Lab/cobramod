#!/usr/bin/env python3
"""Unit-test for package parsing

This test handles the retrieval of data from a local directory or from
the different databases. The module is separated according to the databases:

- TestBiocyc: For Biocyc-related databases
- TestBigg: For BiGG database
- TestKegg: For Kegg
- TestPlantCyc: For PlantCyc
- TestSolcyc: For SolCyc
"""

import unittest
from pathlib import Path

import cobra.core as cobra_core
import cobramod.error as cmod_error
import cobramod.retrieval as cmod_retrieval
import requests
from cobra import __version__ as cobra_version
from cobramod import __version__ as cmod_version
from cobramod.debug import change_to_debug
from cobramod.parsing import bigg as bi
from cobramod.parsing.db_version import DataVersionConfigurator

dir_data = Path(__file__).resolve().parent.joinpath("data")

change_to_debug()
data_conf = DataVersionConfigurator()


# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class TestKegg(unittest.TestCase):
    """
    Parsing and retrival for Kegg
    """

    @classmethod
    def setUpClass(cls):
        data_conf.force_same_version = True

    def test_get_data(self):
        # CASE: Wrong identifier
        self.assertRaises(
            requests.HTTPError,
            cmod_retrieval.get_data,
            directory=dir_data,
            identifier="R08618-wrong",
            database="KEGG",
        )

        # CASE: Regular reaction
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="R08618", database="KEGG"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: Regular compound
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="C01290", database="KEGG"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: RNA
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="GLT-tRNAs", database="pmn:PLANT"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: Module
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="M00001", database="KEGG"
        )
        test_graph = test_data.attributes["pathway"]
        self.assertEqual(first=test_graph["R00200"], second=None)
        self.assertEqual(first=len(test_graph), second=15)

        # CASE: Simple branch
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="M00040", database="KEGG"
        )
        test_graph = test_data.attributes["pathway"]
        self.assertEqual(first=test_graph["R00732"], second=None)
        self.assertEqual(first=test_graph["R00733"], second=None)
        self.assertEqual(first=len(test_graph), second=4)

    def test_parse(self):
        # CASE: Reaction (same as setup)
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="R08618",
            database="KEGG",
            genome="ath",
        )
        self.assertEqual(test_data.attributes["genes"]["rule"], "AT3G19710")
        self.assertEqual(
            test_data.attributes["genes"]["genes"], {"AT3G19710": ""}
        )

        test_reaction = test_data.parse(cobra_core.Model(""), "c")
        self.assertEqual(test_reaction.id, "R08618_c")
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")

        self.assertEqual(first=len(test_reaction.metabolites), second=4)
        self.assertEqual(len(test_reaction.genes), 1)

        # CASE: Compound
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="C01290", database="KEGG"
        )
        test_metabolite = test_data.parse(cobra_core.Model(""), "c")
        self.assertEqual(test_metabolite.id, "C01290_c")

        if not isinstance(test_metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(test_metabolite.formula, "C31H56NO13R")

        # CASE: EC number (not working for now)
        self.assertRaises(
            AttributeError,
            cmod_retrieval.get_data,
            directory=dir_data,
            identifier="7.1.2.2",
            database="KEGG",
        )

        # CASE: Pathway
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="M00001", database="KEGG"
        )
        self.assertEqual("Pathway", test_data.mode)
        self.assertEqual("M00001", test_data.identifier)
        self.assertEqual(15, len(test_data.attributes["pathway"]))


class TestBiocyc(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data_conf.force_same_version = True

    def test_get_data(self):
        # CASE: ID not found
        self.assertRaises(
            requests.HTTPError,
            cmod_retrieval.get_data,
            directory=dir_data,
            identifier="WATER_fake",
            database="META",
        )
        # CASE: Regular META
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="WATER", database="META"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="ACETALD", database="ARA"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="PWY-5337", database="ARA"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: Regular complex lineal pathway

        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="PWY-1187", database="META"
        )

        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 14)
        self.assertCountEqual(test_graph["RXN-2221"], ("RXN-2222", "RXN-2223"))

        # CASE: Complex cyclic pathway

        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="CALVIN-PWY", database="META"
        )
        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 13)
        self.assertCountEqual(
            test_graph["TRIOSEPISOMERIZATION-RXN"],
            ("SEDOBISALDOL-RXN", "F16ALDOLASE-RXN"),
        )

        # CASE: Single-reaction pathway
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="PWY-7344", database="META"
        )

        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 1)
        self.assertEqual(test_graph["UDPGLUCEPIM-RXN"], None)

        # CASE: Super-Pathway
        self.assertRaises(
            cmod_error.SuperpathwayException,
            cmod_retrieval.get_data,
            identifier="PWY-5052",
            database="META",
            directory=dir_data,
        )

    def test_parse(self):
        # CASE: Compound
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="AMP", database="META"
        )
        test_metabolite = test_data.parse(cobra_core.Model(""), "p")

        self.assertEqual(test_metabolite.id, "AMP_p")

        if not isinstance(test_metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(test_metabolite.formula, "C10H12N5O7P1")

        # CASE: Reaction in ARA
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="ARA",
        )

        self.assertEqual(
            test_data.attributes["genes"]["rule"], "AT5G64300 or AT5G59750"
        )
        self.assertCountEqual(
            test_data.attributes["genes"]["genes"].keys(),
            ["AT5G64300", "AT5G59750"],
        )

        # NOTE: retrieve data before this test
        test_reaction = test_data.parse(cobra_core.Model(""), "c")
        self.assertEqual(test_reaction.id, "GTP_CYCLOHYDRO_II_RXN_c")

        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(len(test_reaction.genes), 2)

        # CASE: Reaction in Meta
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="META",
        )
        test_reaction = test_data.parse(cobra_core.Model(""), "c")
        self.assertEqual(test_reaction.id, "GTP_CYCLOHYDRO_II_RXN_c")
        self.assertEqual(test_data.attributes["genes"]["rule"], "")
        self.assertCountEqual(test_data.attributes["genes"]["genes"].keys(), [])
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(len(test_reaction.genes), 0)

        self.assertEqual(test_reaction.bounds, (0, 1000))

        # CASE: Reaction in small subdatabase
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="GCF_000010885",
        )
        self.assertCountEqual(
            test_data.attributes["genes"]["genes"].keys(),
            ["APA22_RS09400"],
        )
        test_reaction = test_data.parse(cobra_core.Model(""), "c")
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(test_reaction.bounds, (0, 1000))

        # CASE: Protein
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="Reduced-hemoproteins",
            database="ARA",
        )

        test_metabolite = test_data.parse(cobra_core.Model(""), "p")

        self.assertEqual(test_metabolite.id, "Reduced_hemoproteins_p")

        if not isinstance(test_metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy Reaction")

        self.assertEqual(test_metabolite.formula, "X")

        # CASE: Pathway
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="PWY-1187",
            database="META",
        )
        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 14)


class TestPlantCyc(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data_conf.force_same_version = True

    def test_get_data(self):
        # CASE: ID not found
        self.assertRaises(
            requests.HTTPError,
            cmod_retrieval.get_data,
            directory=dir_data,
            identifier="WATER_fake",
            database="pmn:PLANT",
        )

        # CASE: Regular PLANT
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="WATER", database="pmn:PLANT"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: META from pmn
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="CPD-15323", database="pmn:META"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: Pathway from CORN
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="CALVIN-PWY", database="pmn:CORN"
        )
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: Regular complex lineal pathway
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="PWY-1187", database="pmn:META"
        )
        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 14)
        self.assertCountEqual(test_graph["RXN-2221"], ("RXN-2222", "RXN-2223"))

        # CASE: Complex cyclic pathway
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="CALVIN-PWY", database="pmn:CORN"
        )
        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 13)
        self.assertCountEqual(
            test_graph["TRIOSEPISOMERIZATION-RXN"],
            ("SEDOBISALDOL-RXN", "F16ALDOLASE-RXN"),
        )

        # CASE: Single-reaction pathway
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="PWY-7344", database="pmn:META"
        )

        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 1)
        self.assertEqual(test_graph["UDPGLUCEPIM-RXN"], None)

        # CASE: Super-Pathway
        self.assertRaises(
            cmod_error.SuperpathwayException,
            cmod_retrieval.get_data,
            identifier="PWY-5052",
            database="pmn:META",
            directory=dir_data,
        )

    def test_parse(self):
        # CASE: Compound
        test_data = cmod_retrieval.get_data(
            directory=dir_data, identifier="AMP", database="pmn:META"
        )
        test_metabolite = test_data.parse(cobra_core.Model(""), "p")

        self.assertEqual(test_metabolite.id, "AMP_p")

        if not isinstance(test_metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(test_metabolite.formula, "C10H12N5O7P1")

        # CASE: Reaction in ARA
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="pmn:ARA",
        )

        self.assertEqual(
            test_data.attributes["genes"]["rule"], "AT5G64300 or AT5G59750"
        )
        self.assertCountEqual(
            test_data.attributes["genes"]["genes"].keys(),
            ["AT5G64300", "AT5G59750"],
        )
        test_reaction = test_data.parse(cobra_core.Model(""), "c")
        self.assertEqual(test_reaction.id, "GTP_CYCLOHYDRO_II_RXN_c")

        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(len(test_reaction.genes), 2)

        # CASE: Reaction in Meta
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="pmn:META",
        )
        test_reaction = test_data.parse(cobra_core.Model(""), "c")
        self.assertEqual(test_reaction.id, "GTP_CYCLOHYDRO_II_RXN_c")
        self.assertEqual(test_data.attributes["genes"]["rule"], "")
        self.assertCountEqual(test_data.attributes["genes"]["genes"].keys(), [])
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(len(test_reaction.genes), 0)

        self.assertEqual(test_reaction.bounds, (0, 1000))

        # CASE: Reaction in another database
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="pmn:SOY",
        )
        test_reaction = test_data.parse(cobra_core.Model(""), "c")
        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(len(test_reaction.genes), 8)

        # CASE: Protein
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="Reduced-hemoproteins",
            database="pmn:ARA",
        )
        test_metabolite = test_data.parse(cobra_core.Model(""), "p")

        self.assertEqual(test_metabolite.id, "Reduced_hemoproteins_p")

        if not isinstance(test_metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy Reaction")

        self.assertEqual(test_metabolite.formula, "X")

        # CASE: Pathway
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="PWY-1187",
            database="pmn:META",
        )
        test_graph = test_data.attributes["pathway"]
        self.assertEqual(len(test_graph), 14)


class TestBigg(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data_conf.force_same_version = True

    def test_find_url(self):
        # CASE: Error
        self.assertRaises(
            requests.HTTPError,
            bi.find_url,
            model_id="e_coli_core",
            query="Invalid",
        )

        # CASE: from specific model
        response, _ = bi.find_url(query="ACALD", model_id="e_coli_core")
        self.assertIsInstance(response, requests.Response)

    def test_get_data(self):
        # CASE: metabolite
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="accoa_c",
            model_id="e_coli_core",
            database="BIGG",
        )
        self.assertEqual(test_data.mode, "Metabolite")
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

        # CASE: Regular reaction from ecoli
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="ACALD",
            model_id="e_coli_core",
            database="BIGG",
        )
        self.assertEqual(test_data.mode, "Reaction")
        self.assertIsInstance(obj=test_data, cls=cmod_retrieval.Data)

    def test_parse(self):
        # CASE: Regular universal reaction
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="ACALDt",
            model_id="universal",
            database="BIGG",
        )
        test_reaction = test_data.parse(cobra_core.Model(""), "c")

        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")
        self.assertEqual(len(test_reaction.genes), 0)

        self.assertEqual(first=len(test_reaction.metabolites), second=2)

        # CASE: Regular universal metabolite
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="coa",
            model_id="universal",
            database="BIGG",
        )
        test_metabolite = test_data.parse(cobra_core.Model(""), "p")
        if not isinstance(test_metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy Reaction")

        self.assertEqual(test_metabolite.formula, "C21H32N7O16P3S")

        # CASE: Ecoli reaction
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="PDH",
            model_id="e_coli_core",
            database="BIGG",
        )
        self.assertEqual(
            test_data.attributes["genes"]["rule"], "b0114 and b0115 and b0116"
        )

        test_reaction = test_data.parse(cobra_core.Model(""), "p")

        if not isinstance(test_reaction, cobra_core.Reaction):
            raise TypeError("Not a valid COBRApy Reaction")

        self.assertEqual(len(test_reaction.genes), 3)

        # CASE: Ecoli metabolite
        test_data = cmod_retrieval.get_data(
            directory=dir_data,
            identifier="co2_c",
            model_id="e_coli_core",
            database="BIGG",
        )
        test_metabolite = test_data.parse(cobra_core.Model(""), "p")
        if not isinstance(test_metabolite, cobra_core.Metabolite):
            raise TypeError("Not a valid COBRApy Reaction")

        self.assertEqual(test_metabolite.id, "co2_p")
        self.assertEqual(test_metabolite.formula, "CO2")


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
