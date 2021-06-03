#!/usr/bin/env python3
"""Unit-test for package parsing

This test handles the retrieval of data from a local directory or from
the different databases. The module is separated according to the databases:

- TestBiocyc: For Biocyc-related databases
- TestBigg: For BiGG database
- TestKegg: For Kegg
"""
from xml.etree.ElementTree import Element
from logging import DEBUG
from pathlib import Path
from unittest import TestCase, main

from requests import HTTPError, Response

from cobramod.debug import debug_log
from cobramod.error import WrongParserError, SuperpathwayWarning
from cobramod.parsing import kegg as kg
from cobramod.parsing import biocyc as bc
from cobramod.parsing import bigg as bi

# Debug must be set in level DEBUG for the test
debug_log.setLevel(DEBUG)
# Setting directory for data
dir_data = Path(__file__).resolve().parent.joinpath("data")
# If data is missing, then do not test. Data should always be the same
if not dir_data.exists():
    raise NotADirectoryError("Data for the test is missing")


class TestKegg(TestCase):
    """
    Parsing and retrival for Kegg
    """

    def test_retrieve_data(self):
        # CASE 0b: Wrong identifier
        self.assertRaises(
            HTTPError,
            kg.retrieve_data,
            directory=dir_data,
            identifier="R08618-wrong",
        )
        #
        # CASE 1: Regular reaction
        test_data = kg.retrieve_data(directory=dir_data, identifier="R08618")
        self.assertIsInstance(obj=test_data, cls=dict)
        # CASE 1: Regular compound
        test_data = kg.retrieve_data(directory=dir_data, identifier="C01290")
        self.assertIsInstance(obj=test_data, cls=dict)

    def test_get_graph(self):
        # CASE 0: Normal Module
        test_dict = kg.retrieve_data(directory=dir_data, identifier="M00001")
        test_graph = kg.get_graph(kegg_dict=test_dict)
        self.assertEqual(first=test_graph["R00200"], second=None)
        self.assertEqual(first=len(test_graph), second=15)
        # CASE 1: Simple branch
        test_dict = kg.retrieve_data(directory=dir_data, identifier="M00040")
        test_graph = kg.get_graph(kegg_dict=test_dict)
        self.assertEqual(first=test_graph["R00732"], second=None)
        self.assertEqual(first=test_graph["R00733"], second=None)
        self.assertEqual(first=len(test_graph), second=4)

    def test__parse_kegg(self):
        # CASE 1a: Reaction (same as setup)
        test_string = kg.retrieve_data(directory=dir_data, identifier="R08618")
        test_dict = kg.KeggParser._parse(
            root=test_string, directory=dir_data, genome="ath"
        )
        self.assertEqual(first=len(test_dict["EQUATION"]), second=4)
        self.assertIsInstance(obj=test_dict["NAME"], cls=str)
        self.assertEqual(first="Reaction", second=test_dict["TYPE"])
        self.assertCountEqual(
            first=test_dict["GENES"]["genes"].keys(), second=["AT3G19710"]
        )
        self.assertEqual(first=test_dict["GENES"]["rule"], second="AT3G19710")
        # CASE 2: Compound
        self.test_string = kg.retrieve_data(
            directory=dir_data, identifier="C01290"
        )
        test_dict = kg.KeggParser._parse(
            root=self.test_string, directory=dir_data
        )

        self.assertEqual(first="Compound", second=test_dict["TYPE"])
        # CASE 2: EC number (not working for now)
        self.test_string = kg.retrieve_data(
            directory=dir_data, identifier="7.1.2.2"
        )
        self.assertRaises(
            NotImplementedError,
            kg.KeggParser._parse,
            root=self.test_string,
            directory=dir_data,
        )
        # CASE 3: Pathway
        test_data = kg.retrieve_data(directory=dir_data, identifier="M00001")
        test_dict = kg._p_pathway(kegg_dict=test_data)
        self.assertEqual(first="Pathway", second=test_dict["TYPE"])
        self.assertEqual(first="M00001", second=test_dict["ENTRY"])


class TestBiocyc(TestCase):
    def test_retrieve_data(self):
        # CASE 1: Directory does not exist
        self.assertRaises(
            NotADirectoryError,
            bc.retrieve_data,
            # args
            directory=Path.cwd().joinpath("noDIr"),
            identifier="WATER",
            database="META",
        )
        # CASE 2: ID not found
        self.assertRaises(
            HTTPError,
            bc.retrieve_data,
            directory=dir_data,
            identifier="WATER_fake",
            database="META",
        )
        # CASE 3: Regular META
        test_element = bc.retrieve_data(
            directory=dir_data, identifier="WATER", database="META"
        )
        self.assertIsInstance(obj=test_element, cls=Element)
        # CASE 4: META from ARA if not available
        test_element = bc.retrieve_data(
            directory=dir_data, identifier="CPD-15323", database="ARA"
        )
        self.assertIsInstance(obj=test_element, cls=Element)

    def test_get_graph(self):
        # CASE 1: Regular complex lineal pathway
        test_root = bc.retrieve_data(
            directory=dir_data, identifier="PWY-1187", database="META"
        )
        test_dict = bc.get_graph(root=test_root)
        self.assertEqual(first=len(test_dict), second=14)
        self.assertCountEqual(
            first=test_dict["RXN-2221"], second=("RXN-2222", "RXN-2223")
        )
        # CASE 2: Complex cyclic pathway
        test_root = bc.retrieve_data(
            directory=dir_data, identifier="CALVIN-PWY", database="META"
        )
        test_dict = bc.get_graph(root=test_root)
        self.assertEqual(first=len(test_dict), second=13)
        self.assertCountEqual(
            first=test_dict["TRIOSEPISOMERIZATION-RXN"],
            second=("SEDOBISALDOL-RXN", "F16ALDOLASE-RXN"),
        )
        # CASE 3: Single-reaction pathway
        test_root = bc.retrieve_data(
            directory=dir_data, identifier="PWY-7344", database="META"
        )
        test_dict = bc.get_graph(root=test_root)
        self.assertEqual(first=len(test_dict), second=1)
        self.assertEqual(first=test_dict["UDPGLUCEPIM-RXN"], second=None)
        # CASE 4: Super-Pathway
        test_root = bc.retrieve_data(
            directory=dir_data, identifier="PWY-5052", database="META"
        )
        self.assertWarns(SuperpathwayWarning, bc.get_graph, root=test_root)

    def test__parse_biocyc(self):
        # CASE 1: Compound
        test_root = bc.retrieve_data(
            directory=dir_data, identifier="AMP", database="META"
        )
        test_dict = bc.BiocycParser._parse(root=test_root, directory=dir_data)
        self.assertEqual(first=test_dict["FORMULA"], second="C10H12N5O7P1")
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 2a: Reaction in ARA
        test_root = bc.retrieve_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="ARA",
        )
        test_dict = bc.BiocycParser._parse(root=test_root, directory=dir_data)
        self.assertEqual(first=len(test_dict["EQUATION"]), second=6)
        self.assertEqual(first=test_dict["EQUATION"]["l_WATER"], second=-3)
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        self.assertCountEqual(
            first=test_dict["GENES"]["genes"].keys(),
            second=["AT5G64300", "AT5G59750"],
        )
        # CASE 2b: Reaction in Meta
        test_root = bc.retrieve_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="META",
        )
        test_dict = bc.BiocycParser._parse(root=test_root, directory=dir_data)
        self.assertCountEqual(
            first=test_dict["GENES"]["genes"].keys(), second=[]
        )
        self.assertEqual(first=test_dict["BOUNDS"], second=(0, 1000))
        # CASE 2c: Reaction in small subdatabase
        test_root = bc.retrieve_data(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="GCF_000010885",
        )
        test_dict = bc.BiocycParser._parse(root=test_root, directory=dir_data)
        self.assertCountEqual(
            first=test_dict["GENES"]["genes"].keys(), second=["APA22_RS09400"]
        )
        self.assertEqual(first=test_dict["BOUNDS"], second=(0, 1000))
        # CASE 3: Protein
        test_root = bc.retrieve_data(
            directory=dir_data,
            identifier="Reduced-hemoproteins",
            database="ARA",
        )
        test_dict = bc.BiocycParser._parse(root=test_root, directory=dir_data)
        self.assertEqual(first=test_dict["TYPE"], second="Protein")
        self.assertEqual(first=test_dict["FORMULA"], second="X")
        # CASE 4: Pathway
        test_root = bc.retrieve_data(
            directory=dir_data, identifier="PWY-1187", database="META"
        )
        test_dict = bc.BiocycParser._parse(root=test_root, directory=dir_data)
        self.assertEqual(first=test_dict["TYPE"], second="Pathway")
        self.assertEqual(first=len(test_dict["PATHWAY"]), second=14)


class TestBigg(TestCase):
    def test__find_url(self):
        # CASE 1: Positive request
        test_response = bi._find_url(
            model_id="e_coli_core", identifier="ACALD"
        )
        self.assertIsInstance(obj=test_response, cls=Response)
        # CASE 2: Negative request
        self.assertRaises(
            HTTPError,
            bi._find_url,
            model_id="e_coli_core",
            identifier="Invalid",
        )

    def test_retrieve_data(self):
        # CASE 0
        test_data = bi.retrieve_data(
            directory=dir_data, identifier="accoa_c", model_id="e_coli_core"
        )
        # CASE 1: Regular reaction from ecoli
        test_data = bi.retrieve_data(
            directory=dir_data, identifier="ACALD", model_id="e_coli_core"
        )
        bi._p_reaction(json_data=test_data)
        self.assertIsInstance(obj=test_data, cls=dict)
        self.assertEqual(
            first=test_data["results"][0]["lower_bound"], second=-1000
        )
        self.assertEqual(first=len(test_data["metabolites"]), second=6)
        # CASE 2: Regular metabolite from universal
        test_data = bi.retrieve_data(
            directory=dir_data, identifier="accoa", model_id="universal"
        )
        self.assertIsInstance(obj=test_data, cls=dict)
        self.assertEqual(
            first=test_data["formulae"][0], second="C23H34N7O17P3S"
        )

    def test__parse_bigg(self):
        # CASE 0: Wrong type
        self.assertRaises(WrongParserError, bi.BiggParser._parse, str())
        # CASE 1: Regular universal reaction
        test_json = bi.retrieve_data(
            directory=dir_data, identifier="ACALD", model_id="universal"
        )
        test_dict = bi.BiggParser._parse(root=test_json)
        self.assertEqual(first=len(test_dict["EQUATION"]), second=6)
        self.assertEqual(first="Reaction", second=test_dict["TYPE"])
        self.assertEqual(
            first={"genes": {}, "rule": ""}, second=test_dict["GENES"]
        )
        # CASE 2: Regular universal metabolite
        test_json = bi.retrieve_data(
            directory=dir_data, identifier="coa", model_id="universal"
        )
        test_dict = bi.BiggParser._parse(root=test_json)
        self.assertEqual(first="Compound", second=test_dict["TYPE"])
        self.assertEqual(first="C21H32N7O16P3S", second=test_dict["FORMULA"])
        # CASE 3: Ecoli reaction
        test_json = bi.retrieve_data(
            directory=dir_data, identifier="PDH", model_id="e_coli_core"
        )
        test_dict = bi.BiggParser._parse(root=test_json)
        self.assertEqual(first=len(test_dict["EQUATION"]), second=6)
        self.assertEqual(first="Reaction", second=test_dict["TYPE"])
        self.assertEqual(
            first="b0114 and b0115 and b0116",
            second=test_dict["GENES"]["rule"],
        )
        self.assertEqual(first=3, second=len(test_dict["GENES"]["genes"]))
        # CASE 4: Ecoli metabolite
        test_json = bi.retrieve_data(
            directory=dir_data, identifier="co2_c", model_id="e_coli_core"
        )
        test_dict = bi.BiggParser._parse(root=test_json)
        self.assertEqual(first="Compound", second=test_dict["TYPE"])
        self.assertEqual(first="CO2", second=test_dict["FORMULA"])


if __name__ == "__main__":
    main(verbosity=2)
