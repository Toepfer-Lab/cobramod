import gzip
import tempfile
from pathlib import Path
from unittest import TestCase
from unittest.mock import MagicMock, patch

import pandas as pd

from cobramod.dbwalker.chebi import Chebi
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable, Uncertain


class TestChebi(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chebi = Chebi()

    def test_name(self):
        self.assertEqual("Chebi", self.chebi.name)

    def test_getStructureFile(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mock_settings = MagicMock()
            mock_settings.cacheDir = Path(tmpdir)

            with patch.object(Chebi, 'settings', mock_settings):
                with patch('cobramod.dbwalker.chebi.requests.get') as mock_get:
                    mock_response = MagicMock()
                    mock_response.status_code = 200

                    tsv_content = "compound_id\tsmiles\tstandard_inchi\tstandard_inchi_key\n"
                    tsv_content += "12345\tCC(=O)O\tInChI=1S/C2H4O2/c1-2(3)4/h4H,1H3\tZKHQWZAMYRWXGA-KQYNXXCUSA-N\n"
                    mock_response.content = gzip.compress(tsv_content.encode())
                    mock_get.return_value = mock_response

                    chebi = Chebi()

                    cache_file = Path(tmpdir) / "chebi" / "chebi-structure.tsv"
                    self.assertTrue(cache_file.exists())
                    self.assertIsInstance(chebi.structure_file, pd.DataFrame)
                    mock_get.assert_called_once()

        # load existing file
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir) / "chebi"
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_file = cache_dir / "chebi-structure.tsv"

            with gzip.open(cache_file, 'wt') as f:
                f.write("compound_id\tsmiles\tstandard_inchi\tstandard_inchi_key\n")
                f.write("12345\tCC(=O)O\tInChI=1S/C2H4O2/c1-2(3)4/h4H,1H3\tZKHQWZAMYRWXGA-KQYNXXCUSA-N\n")

            mock_settings = MagicMock()
            mock_settings.cacheDir = Path(tmpdir)

            with patch.object(Chebi, 'settings', mock_settings):
                with patch('cobramod.dbwalker.chebi.requests.get') as mock_get:
                    chebi = Chebi()

                    mock_get.assert_not_called()
                    self.assertIsNotNone(chebi.structure_file)
                    self.assertEqual(len(chebi.structure_file), 1)
                    self.assertEqual(chebi.structure_file.iloc[0]["compound_id"], 12345)
                    self.assertEqual(chebi.structure_file.iloc[0]["smiles"], "CC(=O)O")
                    self.assertEqual(chebi.structure_file.iloc[0]["standard_inchi"], "InChI=1S/C2H4O2/c1-2(3)4/h4H,1H3")
                    self.assertEqual(chebi.structure_file.iloc[0]["standard_inchi_key"], "ZKHQWZAMYRWXGA-KQYNXXCUSA-N")

    def test_AnnotationPrefix(self):
        self.assertEqual("CHEBI", self.chebi.AnnotationPrefix)

    def test_getGenerellIdentifier(self):
        # 1. Valid ID — returns GenerellIdentifiers with all fields
        result = self.chebi.getGenerellIdentifier("15422")  # known compound
        self.assertIsInstance(result, GenerellIdentifiers)
        self.assertEqual(result.smiles, "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O")
        self.assertEqual(result.inchi, "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1")
        self.assertEqual(result.inchi_key, "ZKHQWZAMYRWXGA-KQYNXXCUSA-N")

        # should return identical entry to the string variant
        result_2 = self.chebi.getGenerellIdentifier(15422)
        self.assertEqual(result, result_2)

        # 2. Non-existent ID — returns Unavailable
        result = self.chebi.getGenerellIdentifier("-999")
        self.assertEqual(result, Unavailable)

        # 3. NaN fields → individual fields become Unavailable
        result = self.chebi.getGenerellIdentifier("37906")
        self.assertIsInstance(result, GenerellIdentifiers)
        self.assertEqual(result.smiles, "*C(=O)C(N)Cc1cncn1")
        self.assertEqual(result.inchi, Unavailable)
        self.assertEqual(result.inchi_key, Unavailable)

        result = self.chebi.getGenerellIdentifier("46622")
        self.assertIsInstance(result, GenerellIdentifiers)
        self.assertEqual(result.smiles, Unavailable)
        self.assertEqual(result.inchi, Unavailable)
        self.assertEqual(result.inchi_key, Unavailable)


        # 4. Duplicate compound_id → raises ValueError
        original_df = self.chebi.structure_file
        dup_df = pd.DataFrame({
            "compound_id": [11111, 11111],
            "smiles": ["C", "CC"],
            "standard_inchi": ["x", "y"],
            "standard_inchi_key": ["a", "b"],
        })
        self.chebi.structure_file = dup_df
        with self.assertRaises(ValueError):
            self.chebi.getGenerellIdentifier("11111")
        self.chebi.structure_file = original_df

    def test_getDBIdentifierFromSmiles(self):
        # 1. Single match — returns compound_id as string
        chebi_id = self.chebi.getDBIdentifierFromSmiles(
            "CC[C@H](C)[C@H](NC(=O)N[C@@H]1CCCCNC(=O)[C@H](COC(C)=O)NC(=O)[C@H](CCc2ccc(O)cc2)N(C)C(=O)[C@H](CCc2ccccc2)NC(=O)[C@H](CCS(C)=O)NC1=O)C(=O)O"  # existing long SMILES
        )
        self.assertEqual("204065", chebi_id)

        # 2. No match — returns Unavailable
        result = self.chebi.getDBIdentifierFromSmiles("XXXXXXXXXXX_NONEXISTENT")
        self.assertEqual(result, Unavailable)

        # 3. Multiple matches — returns Uncertain
        result = self.chebi.getDBIdentifierFromSmiles("*C(=O)O[C@H](CO)COP(=O)([O-])[O-]")
        self.assertIsInstance(result, Uncertain)

        # 4. GenerellIdentifiers input
        gid = GenerellIdentifiers(
            smiles="CCC1CCCC2(CC[C@@H](C)O2)O1",  # same SMILES as test 1
            inchi=Unavailable,
            inchi_key=Unavailable,
        )
        result = self.chebi.getDBIdentifierFromSmiles(gid)
        self.assertEqual("196474", result)

    def test_getDBIdentifierFromInchi(self):
        chebi_id = self.chebi.getDBIdentifierFromInchi(
            "InChI=1S/C15H12O4/c1-11(16)18-14-10-6-5-9-13(14)15(17)19-12-7-3-2-4-8-12/h2-10H,1H3"
        )
        self.assertEqual("135050", chebi_id)

    def test_getDBIdentifierFromInchiKey(self):
        chebi_id = self.chebi.getDBIdentifierFromInchiKey(
            "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"
        )
        self.assertEqual("15422", chebi_id)

    def test_CheckIfCurrentCHEBIFlatFileHasDuplicateCompoundIDs(self):
        duplicates = self.chebi.structure_file[
            self.chebi.structure_file["compound_id"].duplicated(keep=False)
        ]
        self.assertEqual(
            len(duplicates),
            0,
            f"Found {len(duplicates)} rows with duplicate compound_id values: "
            f"{duplicates['compound_id'].unique().tolist()}",
        )
