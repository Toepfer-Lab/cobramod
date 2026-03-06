import json
import tempfile
from pathlib import Path
from unittest import TestCase

from cobramod.dbwalker.cache import Cache, MissmatchError
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable


class TestCache(TestCase):

    def test_addSmiles(self):
        cache = Cache()

        cache.addSmiles("[H]O[H]", "CHEBI:15377")

        self.assertEqual(cache.getBySmiles("[H]O[H]"), set(["CHEBI:15377"]))

        expected = GenerellIdentifiers(smiles="[H]O[H]")
        self.assertEqual(cache.getByID("CHEBI:15377"),expected)
        self.assertEqual(len(cache.smiles_dict), 1)
        self.assertEqual(len(cache.id_dict), 1)

        self.assertEqual(len(cache.inchi_dict), 0)
        self.assertEqual(len(cache.inchi_key_dict), 0)
        self.assertEqual(len(cache._cache_smiles_not_found), 0)
        self.assertEqual(len(cache._cache_inchi_not_found), 0)
        self.assertEqual(len(cache._cache_inchikey_not_found), 0)

        cache.addSmiles("[H]O[H]", Unavailable)

        hit = cache.getBySmiles("[H]O[H]")
        self.assertEqual(hit, Unavailable)

        with self.assertRaises(ValueError):
            cache.addSmiles(123, "CHEBI:15380")

        with self.assertRaises(ValueError):
            cache.addSmiles(None, "CHEBI:15381")

        cache.addSmiles("[H]O[H]", "CHEBI:15377")
        self.assertEqual(cache.getByID("CHEBI:15377"), expected)

        with self.assertRaises(MissmatchError):
            cache.addSmiles("OH", "CHEBI:15377")

        cache_id_dict_size = len(cache.id_dict)
        cache_smiles_dict_size = len(cache.smiles_dict)
        cache_inchi_dict_size = len(cache.inchi_dict)
        cache_inchikey_dict_size = len(cache.inchi_key_dict)
        cache_smiles_not_dict_size = len(cache._cache_smiles_not_found)
        cache_inchi_not_dict_size = len(cache._cache_inchi_not_found)
        cache_inchikey_not_dict_size = len(cache._cache_inchikey_not_found)

        cache.addSmiles("OHCCO", Unavailable)

        self.assertEqual(cache_id_dict_size, len(cache.id_dict))
        self.assertEqual(cache_smiles_dict_size, len(cache.smiles_dict))
        self.assertEqual(cache_inchi_dict_size, len(cache.inchi_dict))
        self.assertEqual(cache_inchikey_dict_size, len(cache.inchi_key_dict))
        self.assertEqual(cache_smiles_not_dict_size + 1, len(cache._cache_smiles_not_found))
        self.assertEqual(cache_inchi_not_dict_size, len(cache._cache_inchi_not_found))
        self.assertEqual(cache_inchikey_not_dict_size, len(cache._cache_inchikey_not_found))

        self.assertEqual(cache.getBySmiles("OHCCO"), Unavailable)

    def test_addInchi(self):
        cache = Cache()

        cache.addInchi("InChI=1S/H2O/h1H2", "CHEBI:15377")

        self.assertEqual(cache.getByInchi("InChI=1S/H2O/h1H2"), set(["CHEBI:15377"]))

        expected = GenerellIdentifiers(inchi="InChI=1S/H2O/h1H2")
        self.assertEqual(cache.getByID("CHEBI:15377"),expected)

        with self.assertRaises(MissmatchError):
            cache.addInchi("InChI=1S/H7O/h1H2", "CHEBI:15377")

        cache.addInchi("InChI=1S/H7O/h1H2", "CHEBI:15378")
        response = cache.getByInchi("InChI=1S/H7O/h1H2")

        self.assertEqual({"CHEBI:15377", "CHEBI:15378"}, response)

    def test_addInchiKey(self):
        cache = Cache()

        cache.addInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "CHEBI:15377")

        self.assertEqual(cache.getByInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N"), set(["CHEBI:15377"]))

        expected = GenerellIdentifiers(inchi_key="XLYOFNOQVPJJNP-UHFFFAOYSA-N")
        self.assertEqual(cache.getByID("CHEBI:15377"),expected)

    def test_addGenerellIdentifiers(self):
        cache = Cache()

        generellIdentifiers = GenerellIdentifiers(
            smiles="[H]O[H]",
            inchi="InChI=1S/H2O/h1H2",
            inchi_key="XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        )
        cache.addGenerellIdentifiers(generellIdentifiers,"CHEBI:15377")

        self.assertEqual(cache.getByID("CHEBI:15377"), generellIdentifiers)


    def test_save_cache(self):

        cache = Cache()
        cache.addSmiles("[H]O[H]", "CHEBI:15377")
        cache.addInchi("InChI=1S/H2O/h1H2", "CHEBI:15377")
        cache.addInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "CHEBI:15377")
        cache.addSmiles("OHCCO", Unavailable)
        cache.addInchi("InChI=1S/MISSING", Unavailable)
        cache.addInchiKey("MISSINGKEY", Unavailable)

        with tempfile.TemporaryDirectory() as tmpdirname:
            cache._cache_folder = Path(tmpdirname)
            cache.save_cache()

            cache_file = Path(tmpdirname) / "cache.xml"
            self.assertTrue(cache_file.exists())

            smiles_not_found_file = Path(tmpdirname) / "smilesNotFound.txt"
            self.assertTrue(smiles_not_found_file.exists())

            inchi_not_found_file = Path(tmpdirname) / "inchiNotFound.txt"
            self.assertTrue(inchi_not_found_file.exists())

            inchikey_not_found_file = Path(tmpdirname) / "inchikeyNotFound.txt"
            self.assertTrue(inchikey_not_found_file.exists())

            with open(smiles_not_found_file, "r") as f:
                smiles_not_found = {line.strip() for line in f if line.strip()}
            self.assertIn("OHCCO", smiles_not_found)

            with open(inchi_not_found_file, "r") as f:
                inchi_not_found = {line.strip() for line in f if line.strip()}
            self.assertIn("InChI=1S/MISSING", inchi_not_found)

            with open(inchikey_not_found_file, "r") as f:
                inchikey_not_found = {line.strip() for line in f if line.strip()}
            self.assertIn("MISSINGKEY", inchikey_not_found)

            with open(cache_file, "r") as f:
                data = json.load(f)
            self.assertIn("CHEBI:15377", data)

    def test_load_cache(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache_dir = Path(tmpdirname)

            cache_data = {
                "CHEBI:15377": {
                    "smiles": "[H]O[H]",
                    "inchi": "InChI=1S/H2O/h1H2",
                    "inchi_key": "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
                }
            }
            with open(cache_dir / "cache.xml", "w") as f:
                json.dump(cache_data, f)

            with open(cache_dir / "smilesNotFound.txt", "w") as f:
                f.write("OHCCO\n")

            with open(cache_dir / "inchiNotFound.txt", "w") as f:
                f.write("InChI=1S/MISSING\n")

            with open(cache_dir / "inchikeyNotFound.txt", "w") as f:
                f.write("MISSINGKEY\n")

            cache = Cache()
            cache._cache_folder = cache_dir
            cache.load_cache()

            self.assertIn("CHEBI:15377", cache.id_dict)
            entry = cache.getByID("CHEBI:15377")
            self.assertEqual(entry.smiles, "[H]O[H]")
            self.assertEqual(entry.inchi, "InChI=1S/H2O/h1H2")
            self.assertEqual(entry.inchi_key, "XLYOFNOQVPJJNP-UHFFFAOYSA-N")

            self.assertEqual(cache.getBySmiles("[H]O[H]"), {"CHEBI:15377"})
            self.assertEqual(cache.getByInchi("InChI=1S/H2O/h1H2"), {"CHEBI:15377"})
            self.assertEqual(cache.getByInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N"), {"CHEBI:15377"})

            self.assertEqual(cache.getBySmiles("OHCCO"), Unavailable)
            self.assertEqual(cache.getByInchi("InChI=1S/MISSING"), Unavailable)
            self.assertEqual(cache.getByInchiKey("MISSINGKEY"), Unavailable)

    def test_save_and_load(self):
        cache = Cache()
        cache.addSmiles("[H]O[H]", "CHEBI:15377")
        cache.addInchi("InChI=1S/H2O/h1H2", "CHEBI:15377")
        cache.addInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "CHEBI:15377")
        cache.addSmiles("CCO", "CHEBI:16236")
        cache.addSmiles("OHCCO", Unavailable)
        cache.addInchi("InChI=1S/MISSING", Unavailable)
        cache.addInchiKey("MISSINGKEY", Unavailable)

        with tempfile.TemporaryDirectory() as tmpdirname:
            cache._cache_folder = Path(tmpdirname)
            cache.save_cache()

            cache2 = Cache()
            cache2._cache_folder = Path(tmpdirname)
            cache2.load_cache()

            self.assertEqual(cache.id_dict.keys(), cache2.id_dict.keys())

            for dbID in cache.id_dict:
                self.assertEqual(cache.getByID(dbID), cache2.getByID(dbID))

            self.assertEqual(cache.smiles_dict.keys(), cache2.smiles_dict.keys())
            self.assertEqual(cache.inchi_dict.keys(), cache2.inchi_dict.keys())
            self.assertEqual(cache.inchi_key_dict.keys(), cache2.inchi_key_dict.keys())

            self.assertEqual(cache._cache_smiles_not_found, cache2._cache_smiles_not_found)
            self.assertEqual(cache._cache_inchi_not_found, cache2._cache_inchi_not_found)
            self.assertEqual(cache._cache_inchikey_not_found, cache2._cache_inchikey_not_found)

            self.assertEqual(cache2.getBySmiles("OHCCO"), Unavailable)
            self.assertEqual(cache2.getByInchi("InChI=1S/MISSING"), Unavailable)
            self.assertEqual(cache2.getByInchiKey("MISSINGKEY"), Unavailable)

            self.assertEqual(cache2.getBySmiles("[H]O[H]"), {"CHEBI:15377"})
            self.assertEqual(cache2.getBySmiles("CCO"), {"CHEBI:16236"})

    def test_cache_init_loads_from_dir(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache_dir = Path(tmpdirname)

            cache_data = {
                "CHEBI:15377": {
                    "smiles": "[H]O[H]",
                    "inchi": "InChI=1S/H2O/h1H2",
                    "inchi_key": "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
                }
            }
            with open(cache_dir / "cache.xml", "w") as f:
                json.dump(cache_data, f)

            with open(cache_dir / "smilesNotFound.txt", "w") as f:
                f.write("")

            with open(cache_dir / "inchiNotFound.txt", "w") as f:
                f.write("")

            with open(cache_dir / "inchikeyNotFound.txt", "w") as f:
                f.write("")

            cache = Cache(cache_dir=cache_dir)

            self.assertIn("CHEBI:15377", cache.id_dict)
            entry = cache.getByID("CHEBI:15377")
            self.assertEqual(entry.smiles, "[H]O[H]")

    def test_cache_dicts_in_sync(self):
        cache = Cache()
        cache.addSmiles("[H]O[H]", "C:1")
        cache.addInchi("InChI=1S/H2O/h1H2", "C:1")
        cache.addInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "C:1")

        entry = cache.getByID("C:1")
        self.assertEqual(entry.smiles, "[H]O[H]")
        self.assertEqual(entry.inchi, "InChI=1S/H2O/h1H2")
        self.assertEqual(entry.inchi_key, "XLYOFNOQVPJJNP-UHFFFAOYSA-N")

        cache.addInchiKey("DEADBEEF-KEY", Unavailable)
        self.assertIs(cache.getByInchiKey("DEADBEEF-KEY"), Unavailable)
