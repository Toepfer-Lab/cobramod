from unittest import TestCase

from cobramod.dbwalker.pubchem import PubChem


class TestGetPubChem(TestCase):
    def setUp(self):
        self.pubchem = PubChem()

    def test_getGenerellIdentifier(self):
        pubchemID = "12345"

        identifier = self.pubchem.getGenerellIdentifier(pubchemID)
        print(identifier)
