from unittest import TestCase
from unittest.mock import MagicMock, patch

from cobramod.dbwalker.pubchem import PubChem


class TestGetPubChem(TestCase):
    def setUp(self):
        self.pubchem = PubChem()

    def test_name(self):
        self.assertEqual("PubChem", self.pubchem.name)

    def test_AnnotationPrefix(self):
        self.assertEqual("pubchem.compound", self.pubchem.AnnotationPrefix)

   # def test_getGenerellIdentifier(self):

   # def test_getDBIdentifierFromSmiles(self):

   # def test_getDBIdentifierFromInchi(self):

    def test_getDBIdentifierFromInchiKey(self):
        mock_session = MagicMock()

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = '945872\n'
        mock_session.post.return_value = mock_response

        # normal case
        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromInchiKey("NELACKPYIURSBY-MRXNPFEDSA-N")
            self.assertEqual("945872", pubchemID)

        # does not exist
        # multiple results
        mock_response.text = '1264\n102212881\n'
        mock_session.post.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromInchiKey("PFRKDKQPQSBYQX-UHFFFAOYSA-O")
            self.assertEqual("945872", pubchemID)

        # live test
        self.pubchem._cache.clear_cache()
        pubchemID = self.pubchem.getDBIdentifierFromInchiKey("JDSRHVWSAMTSSN-UHFFFAOYSA-N ")
        self.assertEqual("1426", pubchemID)
