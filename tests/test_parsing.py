#!/usr/bin/env python3
from logging import DEBUG
from pathlib import Path
import unittest

from cobramod.debug import debug_log
import cobramod.parsing.kegg as kg
import cobramod.parsing.biocyc as bc
import cobramod.parsing.bigg as bi

debug_log.setLevel(DEBUG)
dir_input = Path.cwd().joinpath("tests").joinpath("input")
dir_data = Path.cwd().joinpath("tests").joinpath("data")


if not dir_data.exists():
    dir_data.mkdir(parents=True)


class ParsingTesting(unittest.TestCase):
    def setUp(self):
        self.raw = """
ENTRY       R08618                      Reaction
NAME        L-methionine:2-oxo-acid aminotransferase
DEFINITION  L-Methionine + 2-Oxo acid <=> 4-Methylthio-2-\
oxobutanoic acid + L-Amino acid
EQUATION    C00073 + C00161 <=> C01180 + C00151
RCLASS      RC00006  C00073_C01180
            RC00025  C00151_C00161
ENZYME      2.6.1.88
PATHWAY     rn00966  Glucosinolate biosynthesis
            rn01110  Biosynthesis of secondary metabolites
            rn01210  2-Oxocarboxylic acid metabolism
ORTHOLOGY   K14287  methionine transaminase [EC:2.6.1.88]
            K21346  methionine transaminase [EC:2.6.1.88]
            ///"""

    def test__get_unformatted_kegg(self):
        # CASE 0a: Colon (:) in identifier (not implemented yet)
        self.assertRaises(
            NotImplementedError,
            kg._get_unformatted_kegg,
            directory=dir_data,
            identifier="rn:R08618",
        )
        # CASE 0b: Wrong identifier
        self.assertRaises(
            Warning,
            kg._get_unformatted_kegg,
            directory=dir_data,
            identifier="R08618.",
        )
        #
        # CASE 1: Regular reaction
        test_data = kg._get_unformatted_kegg(
            directory=dir_data, identifier="R08618"
        )
        self.assertIsInstance(obj=test_data, cls=str)
        # CASE 1: Regular compound
        test_data = kg._get_unformatted_kegg(
            directory=dir_data, identifier="C01290"
        )
        self.assertIsInstance(obj=test_data, cls=str)

    def test__parse_kegg(self):
        # CASE 1: Reaction (same as setUP)
        test_dict = kg.KeggParser._parse(raw=self.raw)
        self.assertEqual(first=len(test_dict["EQUATION"]), second=4)
        self.assertIsInstance(obj=test_dict["NAME"], cls=str)
        self.assertEqual(first="Reaction", second=test_dict["TYPE"])
        # CASE 2: Compound
        self.raw = kg._get_unformatted_kegg(
            directory=dir_data, identifier="C01290"
        )
        test_dict = kg.KeggParser._parse(raw=self.raw)

        self.assertEqual(first="Compound", second=test_dict["TYPE"])
        # CASE 2: EC number (not working for now)
        self.raw = kg._get_unformatted_kegg(
            directory=dir_data, identifier="7.1.2.2"
        )
        self.assertRaises(
            NotImplementedError, kg.KeggParser._parse, raw=self.raw
        )
        # CASE 3: Pathway
        test_raw = kg._get_unformatted_kegg(
            directory=dir_data, identifier="M00001"
        )
        test_data = kg._create_dict(raw=test_raw)
        test_dict = kg._p_pathway(kegg_dict=test_data)
        self.assertEqual(first="Pathway", second=test_dict["TYPE"])
        self.assertEqual(first="M00001", second=test_dict["ENTRY"])

    def test__get_xml_from_biocyc(self):
        # CASE 1: Directory does not exist
        self.assertRaises(
            NotADirectoryError,
            bc._get_xml_from_biocyc,
            # args
            directory=Path.cwd().joinpath("noDIr"),
            identifier="WATER",
            database="META",
        )
        # CASE 2: ID not found
        self.assertRaises(
            Warning,
            bc._get_xml_from_biocyc,
            directory=dir_data,
            identifier="WATE",
            database="META",
        )

    def test__parse_biocyc(self):
        # CASE 1: Compound
        test_root = bc._get_xml_from_biocyc(
            directory=dir_data, identifier="AMP", database="META"
        )
        test_dict = bc.BiocycParser._parse(root=test_root)
        self.assertEqual(first=test_dict["FORMULA"], second="C10H12N5O7P1")
        self.assertEqual(first=test_dict["TYPE"], second="Compound")
        # CASE 2: Reaction
        test_root = bc._get_xml_from_biocyc(
            directory=dir_data,
            identifier="GTP-CYCLOHYDRO-II-RXN",
            database="META",
        )
        test_dict = bc.BiocycParser._parse(root=test_root)
        self.assertEqual(first=len(test_dict["EQUATION"]), second=6)
        self.assertEqual(first=test_dict["EQUATION"]["l_WATER"], second=-3)
        self.assertEqual(first=test_dict["TYPE"], second="Reaction")
        self.assertEqual(first=test_dict["BOUNDS"], second=(0, 1000))
        # CASE 3: Protein
        test_root = bc._get_xml_from_biocyc(
            directory=dir_data,
            identifier="Reduced-hemoproteins",
            database="ARA",
        )
        test_dict = bc.BiocycParser._parse(root=test_root)
        self.assertEqual(first=test_dict["TYPE"], second="Protein")
        self.assertEqual(first=test_dict["FORMULA"], second="X")
        # CASE 4: Pathway
        test_root = bc._get_xml_from_biocyc(
            directory=dir_data, identifier="PWY-1187", database="META"
        )
        test_dict = bc.BiocycParser._parse(root=test_root)
        self.assertEqual(first=test_dict["TYPE"], second="Pathway")
        self.assertEqual(first=len(test_dict["PATHWAY"]), second=13)
        self.assertEqual(first=len(test_dict["SET"]), second=14)

    def test__get_json_bigg(self):
        # CASE 0
        test_data = bi._get_json_bigg(
            directory=dir_data,
            identifier="accoa_c",
            model_id="e_coli_core",
            object_type="metabolite",
        )
        # CASE 1: Regular reaction from ecoli
        test_data = bi._get_json_bigg(
            directory=dir_data,
            identifier="ACALD",
            model_id="e_coli_core",
            object_type="reaction",
        )
        bi._p_reaction(json_data=test_data)
        self.assertIsInstance(obj=test_data, cls=dict)
        self.assertEqual(
            first=test_data["results"][0]["lower_bound"], second=-1000
        )
        self.assertEqual(first=len(test_data["metabolites"]), second=6)
        # CASE 2: Regular metabolite from universal
        test_data = bi._get_json_bigg(
            directory=dir_data,
            identifier="accoa",
            model_id="universal",
            object_type="metabolite",
        )
        self.assertIsInstance(obj=test_data, cls=dict)
        self.assertEqual(
            first=test_data["formulae"][0], second="C23H34N7O17P3S"
        )

    def test__parse_bigg(self):
        # CASE 1: Regular universal reaction
        test_json = bi._get_json_bigg(
            directory=dir_data,
            identifier="ACALD",
            model_id="universal",
            object_type="reaction",
        )
        test_dict = bi.BiggParser._parse(json_data=test_json)
        self.assertEqual(first=len(test_dict["EQUATION"]), second=6)
        self.assertEqual(first="Reaction", second=test_dict["TYPE"])
        # CASE 2: Regular universal metabolite
        test_json = bi._get_json_bigg(
            directory=dir_data,
            identifier="coa",
            model_id="universal",
            object_type="metabolite",
        )
        test_dict = bi.BiggParser._parse(json_data=test_json)
        self.assertEqual(first="Compound", second=test_dict["TYPE"])
        self.assertEqual(first="C21H32N7O16P3S", second=test_dict["FORMULA"])
        # CASE 3: Ecoli reaction
        test_json = bi._get_json_bigg(
            directory=dir_data,
            identifier="PDH",
            model_id="e_coli_core",
            object_type="reaction",
        )
        test_dict = bi.BiggParser._parse(json_data=test_json)
        self.assertEqual(first=len(test_dict["EQUATION"]), second=6)
        self.assertEqual(first="Reaction", second=test_dict["TYPE"])
        # CASE 4: Ecoli metabolite
        test_json = bi._get_json_bigg(
            directory=dir_data,
            identifier="co2_c",
            model_id="e_coli_core",
            object_type="metabolite",
        )
        test_dict = bi.BiggParser._parse(json_data=test_json)
        self.assertEqual(first="Compound", second=test_dict["TYPE"])
        self.assertEqual(first="CO2", second=test_dict["FORMULA"])

    def test__get_db_names(self):
        test_dict = kg._create_dict(
            raw=kg._get_unformatted_kegg(
                directory=dir_data, identifier="C00073"
            )
        )
        test_string = kg._get_db_names(data_dict=test_dict, pattern="ChEBI")
        self.assertEqual(first=test_string, second="16643")


if __name__ == "__main__":
    unittest.main(verbosity=2)
