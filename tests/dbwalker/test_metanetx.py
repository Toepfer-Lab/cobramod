from unittest import TestCase

from cobramod import Settings
from cobramod.dbwalker.dataclasses import (
    GenerellIdentifiers,
    Unavailable,
    Uncertain,
)
from cobramod.dbwalker.metanetx import MetaNetX


class TestMetaNetX(TestCase):
    @classmethod
    def setUpClass(cls):
        # Downloads files on first run, reuses cached files afterwards
        cls.mnx = MetaNetX()

    def setUp(self):
        self.mnx = self.__class__.mnx

    def test_name(self):
        self.assertEqual("MetaNetX", self.mnx.name)

    def test_AnnotationPrefix(self):
        self.assertEqual("metanetx.chemical", self.mnx.AnnotationPrefix)

    # ── getGenerellIdentifier ────────────────────────────────────────

    def test_getGenerellIdentifier(self):
        gid = self.mnx.getGenerellIdentifier("MNXM1051475")

        expected_gid = GenerellIdentifiers(
            smiles="CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCCC(=O)O[C@H](COC(=O)CCCCCCCC/C=C\C/C=C\C/C=C\CCCCC)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=CC(N)=NC2=O)[C@H](O)[C@@H]1O",
            inchi_key="WEYJPVGCJYJVMD-YQSRUJTQSA-L",
            inchi="InChI=1S/C56H91N3O15P2/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-35-37-39-41-51(60)69-45-48(72-52(61)42-40-38-36-34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2)46-70-75(65,66)74-76(67,68)71-47-49-53(62)54(63)55(73-49)59-44-43-50(57)58-56(59)64/h11-14,17-20,23-26,30,32,43-44,48-49,53-55,62-63H,3-10,15-16,21-22,27-29,31,33-42,45-47H2,1-2H3,(H,65,66)(H,67,68)(H2,57,58,64)/p-2/b13-11-,14-12-,19-17-,20-18-,25-23-,26-24-,32-30-/t48-,49-,53-,54-,55-/m1/s1"
        )

        self.assertTrue(expected_gid == gid)

    def test_getGenerellIdentifier_not_found(self):
        result = self.mnx.getGenerellIdentifier("MNXM_DOES_NOT_EXIST_99999999")
        self.assertIs(Unavailable, result)

    # ── getDBIdentifierFromSmiles ────────────────────────────────────

    def test_getDBIdentifierFromSmiles(self):
        result = self.mnx.getDBIdentifierFromSmiles("CC1=C(C)SC(N2[C@@H](C#N)[C@@H](C3=CC=CC=C3)[C@@H]2CO)=N1")
        expected_ID = "MNXM170624"

        self.assertEqual(expected_ID, result)
        

    def test_getDBIdentifierFromSmiles_not_found(self):
        result = self.mnx.getDBIdentifierFromSmiles("THIS_IS_NOT_A_REAL_SMILES")
        self.assertIs(Unavailable, result)

    def test_getDBIdentifierFromSmiles_with_GenerellIdentifiers(self):
        gid = GenerellIdentifiers(smiles="N[C@H](CSSC[C@@H](N)C(=O)O)C(=O)O")
        result = self.mnx.getDBIdentifierFromSmiles(gid)

        self.assertEqual(result, "MNXM17062")

    # ── getDBIdentifierFromInchi ─────────────────────────────────────

    def test_getDBIdentifierFromInchi(self):
        result = self.mnx.getDBIdentifierFromInchi("InChI=1S/C18H15FN2O2/c19-14-9-5-4-8-13(14)18(23)21-15(10-20)17(16(21)11-22)12-6-2-1-3-7-12/h1-9,15-17,22H,11H2/t15-,16+,17+/m0/s1")
        
        self.assertEqual("MNXM170669",result)

    def test_getDBIdentifierFromInchi_not_found(self):
        result = self.mnx.getDBIdentifierFromInchi("InChI=1S/DOES_NOT_EXIST")
        self.assertIs(Unavailable, result)

    def test_getDBIdentifierFromInchi_with_GenerellIdentifiers(self):
        gid = GenerellIdentifiers(inchi="InChI=1S/C24H40O3/c1-15(7-10-21(26)27)17-8-9-18-22-19(11-13-24(17,18)3)23(2)12-5-4-6-16(23)14-20(22)25/h15-20,22,25H,4-14H2,1-3H3,(H,26,27)/t15-,16+,17-,18+,19+,20+,22+,23+,24-/m1/s1")
        result = self.mnx.getDBIdentifierFromInchi(gid)
        self.assertIsNotNone(result)
        self.assertIsNot(Unavailable, result)

    # ── getDBIdentifierFromInchiKey ──────────────────────────────────

    def test_getDBIdentifierFromInchiKey_water(self):
        # Standard InChIKey for water
        result = self.mnx.getDBIdentifierFromInchiKey("PGRXUSLZBWAWNM-HGLLRIGKSA-N")
        self.assertEqual("MNXM244989", result)

    def test_getDBIdentifierFromInchiKey_not_found(self):
        result = self.mnx.getDBIdentifierFromInchiKey("FAKEKEY-DOESNOTEXIST-N")
        self.assertIs(Unavailable, result)

    def test_getDBIdentifierFromInchiKey_with_GenerellIdentifiers(self):
        gid = GenerellIdentifiers(inchi_key="LDSYMZUTAQHRHU-FGLMURGOSA-M")
        result = self.mnx.getDBIdentifierFromInchiKey(gid)
        self.assertEqual( "MNXM383513", result)

    # ── Round-trip consistency ───────────────────────────────────────

    def test_round_trip(self):
        # Get identifiers for water, then look up MNX ID from each
        original_id = "MNXM369457"
        gid = self.mnx.getGenerellIdentifier(original_id)

        ecpected_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)O[C@@H](COC(=O)CCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC"
        self.assertEqual(gid.smiles, ecpected_smiles)

        result = self.mnx.getDBIdentifierFromSmiles(gid.smiles)
        self.assertEqual(original_id, result)

        expected_inchi = "InChI=1S/C54H104O6/c1-4-7-10-13-16-19-22-24-25-26-27-28-29-31-33-36-39-42-45-48-54(57)60-51(49-58-52(55)46-43-40-37-34-21-18-15-12-9-6-3)50-59-53(56)47-44-41-38-35-32-30-23-20-17-14-11-8-5-2/h51H,4-50H2,1-3H3/t51-/m0/s1"
        self.assertEqual(expected_inchi, gid.inchi)

        result = self.mnx.getDBIdentifierFromInchi(gid.inchi)
        self.assertEqual(result, original_id)


        expected_inchikey = "WYDOMVHRQQCHML-XHIZWQFQSA-N"
        self.assertEqual(expected_inchikey, gid.inchi_key)

        result = self.mnx.getDBIdentifierFromInchiKey(gid.inchi_key)
        self.assertEqual(result, original_id)

    def test_save_cache(self):
        # save_cache is a no-op, just verify it doesn't raise
        self.mnx.save_cache()

    # ── Header parsing ───────────────────────────────────────────────

    def test_parse_header_chem_prop(self):
        cache_dir = Settings().cacheDir / "MetaNetX"
        skip, cols = MetaNetX._MetaNetX__parse_header(cache_dir / "chem_prop.tsv")

        expected_cols = ["ID", "name", "reference", "formula", "charge", "mass", "InChI", "InChIKey", "SMILES"]
        self.assertEqual(expected_cols, cols)
        self.assertEqual(skip, 389)
        self.assertGreater(skip, 0, "There should be at least one comment line")

    def test_parse_header_chem_xref(self):
        cache_dir = Settings().cacheDir / "MetaNetX"
        skip, cols = MetaNetX._MetaNetX__parse_header(cache_dir / "chem_xref.tsv")

        expected_cols = ["source", "ID", "description"]
        self.assertEqual(expected_cols, cols)
        self.assertGreater(skip, 0, "There should be at least one comment line")
