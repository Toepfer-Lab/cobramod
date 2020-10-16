#!/usr/bin/env python3
from pathlib import Path
import unittest
import cobramod.parsing.kegg as kg
import cobramod.parsing.biocyc as bc
from cobramod.debug import debug_log
from logging import DEBUG
import xml.etree.ElementTree as ET

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
        # CASE 0: Colon (:) in identifier (not implemented yet)
        self.assertRaises(
            NotImplementedError,
            kg._get_unformatted_kegg,
            directory=dir_data,
            identifier="rn:R08618",
        )
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
        self.assertEqual(first="C01290", second=test_dict["ENTRY"])
        self.assertEqual(first="Compound", second=test_dict["TYPE"])
        # CASE 2: EC number (not working for now)
        self.raw = kg._get_unformatted_kegg(
            directory=dir_data, identifier="7.1.2.2"
        )
        self.assertRaises(
            NotImplementedError, kg.KeggParser._parse, raw=self.raw
        )
        # CASE 3: Pathway (Not implemented)
        self.raw = kg._get_unformatted_kegg(
            directory=dir_data, identifier="ath00966"
        )
        self.assertRaises(
            NotImplementedError, kg.KeggParser._parse, raw=self.raw
        )

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
        # CASE 3: Proper retrieval.
        dir_data.joinpath("META").joinpath("WATER.xml").unlink()
        self.assertIsInstance(
            bc._get_xml_from_biocyc(
                directory=dir_data, identifier="WATER", database="META"
            ),
            ET.Element,
        )
        # CASE 4: not found in database
        if dir_data.joinpath("CPD-15326.xml").exists():
            dir_data.joinpath("CPD-15326.xml").unlink()
        self.assertRaises(
            Warning,
            bc._get_xml_from_biocyc,
            directory=dir_data,
            identifier="CPD-15326",
            database="ARA",
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
        self.assertEqual(first=test_dict["EQUATION"]["WATER"], second=-3)
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
