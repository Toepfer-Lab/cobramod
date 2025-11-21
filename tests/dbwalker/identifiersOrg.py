import logging
from unittest import TestCase
from unittest.mock import patch, Mock

from cobramod.dbwalker.identifiersORG import Validator

logger = logging.getLogger("cobramod")
logger.setLevel(logging.DEBUG)


class TestValidator(TestCase):
    def test_initialisation(self):
        validator = Validator()

        self.assertTrue(isinstance(validator, Validator))

    @patch("requests.get")
    def test_second_initialisation(self, mocked_get):
        # Test that the singleton works and therefore calls to identifiers.org are minimized.

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = {}
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            validator = Validator()


            self.assertTrue(isinstance(validator, Validator))
            self.assertGreater(len(validator.lookup), 0)

    def test_Bigg(self):
        validator = Validator()

        self.assertTrue(validator.validate_id_pattern("bigg.metabolite", "Glc_aD"))

        self.assertFalse(validator.validate_id_pattern("bigg.metabolite", ""))

        self.assertFalse(validator.validate_id_pattern("bigg.metabolite", "-.AT"))

    def test_BioCyc(self):
        validator = Validator()

        self.assertFalse(validator.validate_id_pattern("biocyc", "asdwä"))
        self.assertFalse(validator.validate_id_pattern("biocyc", ""))
        self.assertFalse(validator.validate_id_pattern("biocyc", "ecoli:CPD0-1029"))

        self.assertTrue(validator.validate_id_pattern("biocyc", "ECOLI:CPD0-1029"))
        self.assertTrue(validator.validate_id_pattern("biocyc", "CPD0-1029"))