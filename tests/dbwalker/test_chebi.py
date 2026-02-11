from unittest import TestCase

from cobramod.dbwalker.chebi import Chebi


class TestChebi(TestCase):
    def setUp(self):
        self.chebi = Chebi()

    def test_name(self):
        self.assertEqual("Chebi", self.chebi.name)

    def test_AnnotationPrefix(self):
        self.assertEqual("CHEBI", self.chebi.AnnotationPrefix)

    def test_getDBIdentifierFromSmiles(self):

        chebi_id = self.chebi.getDBIdentifierFromSmiles("CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O[C@@H]5O[C@H](CO)[C@H](O)[C@H](O[C@]6(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)[C@H]([C@H](O)[C@H](O)CO)O6)[C@H]5O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](CO[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)O[C@H]1O")
        self.assertEqual("156845", chebi_id)

    def test_getDBIdentifierFromInchi(self):

        chebi_id = self.chebi.getDBIdentifierFromInchi("InChI=1S/C15H12O4/c1-11(16)18-14-10-6-5-9-13(14)15(17)19-12-7-3-2-4-8-12/h2-10H,1H3")
        self.assertEqual("135050", chebi_id)

    def test_getDBIdentifierFromInchiKey(self):

        chebi_id = self.chebi.getDBIdentifierFromInchiKey("ZKHQWZAMYRWXGA-KQYNXXCUSA-N")
        self.assertEqual("15422", chebi_id)