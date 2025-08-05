import tempfile
from pathlib import Path
from unittest import TestCase


from cobramod.core.crossreferences import validate_id
from cobramod.debug import change_to_debug

change_to_debug()

class TestCrossReferencesPubChem(TestCase):
    def test_validate_id(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)

            XRef_existing = "chebi:CHEBI:27732"
            XREF_non_existing = "chebi:CHEBI:0"
            Wrong_provider = "NonExisting:123456"


            result = validate_id(XRef_existing)

            self.assertTrue(result)


