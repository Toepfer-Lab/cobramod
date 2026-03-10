import gzip
import tempfile
from pathlib import Path
from unittest import TestCase
from unittest.mock import MagicMock, patch

import pandas as pd

from cobramod.dbwalker.chebi import Chebi


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

    def test_getDBIdentifierFromSmiles(self):
        chebi_id = self.chebi.getDBIdentifierFromSmiles(
            "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O[C@@H]5O[C@H](CO)[C@H](O)[C@H](O[C@]6(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)[C@H]([C@H](O)[C@H](O)CO)O6)[C@H]5O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](CO[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)O[C@H]1O"
        )
        self.assertEqual("156845", chebi_id)

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
