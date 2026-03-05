import logging
import shutil
import tempfile
from unittest import TestCase
from unittest.mock import MagicMock, Mock, PropertyMock, patch

import requests
from requests import Session

from cobramod.dbwalker.BioCyc import BioCyc
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable
from cobramod.settings import Settings

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

logger = logging.getLogger("cobramod")
logger.propagate = True
logger.setLevel(logging.DEBUG)


class TestGetCompoundInfoByBiocycId(TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.biocyc = BioCyc(cachedir=self.test_dir)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_getGeneralIdentifiers(self):
        biocycID = "ECOLI:LACTALD"

        gID = self.biocyc.getGenerellIdentifier(biocycID)

        self.assertIsInstance(gID, GenerellIdentifiers)
        self.assertEqual("BSABBBMNWQWLLU-VKHMYHEASA-N", gID.inchi_key)
        self.assertEqual(
            "InChI=1S/C3H6O2/c1-3(5)2-4/h2-3,5H,1H3/t3-/m0/s1", gID.inchi
        )
        self.assertEqual("C[C@@H](C=O)O", gID.smiles)

        gID = self.biocyc.getGenerellIdentifier(biocycID)

        self.assertIsInstance(gID, GenerellIdentifiers)
        self.assertEqual("BSABBBMNWQWLLU-VKHMYHEASA-N", gID.inchi_key)
        self.assertEqual(
            "InChI=1S/C3H6O2/c1-3(5)2-4/h2-3,5H,1H3/t3-/m0/s1", gID.inchi
        )
        self.assertEqual("C[C@@H](C=O)O", gID.smiles)

    def test_get_BioCycIDviaSmiles(self):
        # simple smiles but not unique
        biocycID = self.biocyc.getDBIdentifierFromSmiles("O")
        self.assertEqual(biocycID, Unavailable)

        # identifier not present
        smiles = "CC(=O)OCOC(=O)C"
        biocycID = self.biocyc.getDBIdentifierFromSmiles(smiles)
        self.assertEqual(biocycID, Unavailable)

        # complicated id
        smiles = "C(O)[C@H]2(O[C@@H](O)[C@H](O[C@@H]1([C@H](O)[C@H]([C@H](O)[C@@H](CO)O1)O))[C@@H](O)[C@H](O)2)"
        biocycID = self.biocyc.getDBIdentifierFromSmiles(smiles)
        self.assertEqual(first=biocycID, second="META:BETA-KOJIBIOSE")

        smiles = r"C(OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@H]3(O[C@@H](N1(C2(\C(\N=C/1)=C(N)/N=C\N=2)))[C@H](O)[C@H](O)3)"
        biocycID = self.biocyc.getDBIdentifierFromSmiles(smiles)
        self.assertEqual(first=biocycID, second="META:ATP")

    def test_get_BioCycIDviaInChI(self):
        mock_session = MagicMock()

        mock_response = MagicMock()
        mock_response.status_code = 200

        self.biocyc.clear_caches()
        mock_session.post.return_value.json.return_value = {
            "RESULTS": [
                {
                    "OBJECT-ID": "WATER",
                    "COMMON-NAME": "H<sub>2</sub>O",
                    "INCHI-STRING": "InChI=1S/H2O/h1H2",
                    "INCHI-KEYS": None,
                }
            ]
        }

        with patch.object(
                Settings, "_biocycSession", new_callable=PropertyMock
        ) as mock_prop:
            mock_prop.return_value = mock_session

            biocycID = self.biocyc.getDBIdentifierFromInchi("InChI=1S/H2O/h1H2")
            self.assertEqual("META:WATER", biocycID)

        # Not existing inchi
        self.biocyc.clear_caches()
        mock_session.post.return_value.json.return_value = {"RESULTS": None}

        with patch.object(
                Settings, "_biocycSession", new_callable=PropertyMock
        ) as mock_prop:
            mock_prop.return_value = mock_session

            biocycID = self.biocyc.getDBIdentifierFromInchi(
                "InChI=1S/C3H9O6P/c4-1-3(5)2-9-10(6,7)8/h3-5H,1-2H2,(H2,6,7,8) "
            )
            self.assertEqual(Unavailable, biocycID)

        # Live test
        self.biocyc.clear_caches()
        biocycID = self.biocyc.getDBIdentifierFromInchi("InChI=1S/H2O/h1H2")
        self.assertEqual("META:WATER", biocycID)

    def test_get_BioCycIDviaInChIKey(self):
        mock_session = MagicMock()

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = 'Successful queries\nQuery\tBioCyc Common-Name\tBioCyc\tKegg\tCHEBI\tPubChem\tHMDB\tChemSpider\tMetabolights\tMetanetX\tBiGG\tSEED\tInChi\t\nInChIKey:GUBGYTABKSRVRQ-DCSYEGIMSA-N\tbeta-lactose\t<A NAME="Sbeta-lactose" HREF="/compound?orgid=META&id=CPD-15971">&beta;-lactose</A>\tC01970\t36218\t6134\tHMDB41627\tNIL\tMTBLC36218\tNIL\tNIL\tNIL\n\nAmbigous Queries\nNone\n\nUnknown Queries\nNone\n'
        mock_session.post.return_value = mock_response

        with patch.object(
                Settings, "_biocycSession", new_callable=PropertyMock
        ) as mock_prop:
            mock_prop.return_value = mock_session

            self.biocyc.clear_caches()
            biocycID = self.biocyc.getDBIdentifierFromInchiKey(
                "GUBGYTABKSRVRQ-DCSYEGIMSA-N"
            )
            self.assertEqual("META:CPD-15971", biocycID)

        mock_response.text = 'Successful queries\nQuery\tBioCyc Common-Name\tBioCyc\tKegg\tCHEBI\tPubChem\tHMDB\tChemSpider\tMetabolights\tMetanetX\tBiGG\tSEED\tInChi\t\nInChIKey:GUBGYTABKSRVRQ-ASMJPISFSA-N\talpha-maltose\t<A NAME="Salpha-maltose" HREF="/compound?orgid=META&id=ALPHA-MALTOSE">&alpha;-maltose</A>\tC00897\t18167\t439341\tHMDB00163\t388469\tMTBLC18167\tNIL\tNIL\tNIL\n\nAmbigous Queries\nNone\n\nUnknown Queries\nNone\n'
        mock_session.post.return_value = mock_response

        with patch.object(
                Settings, "_biocycSession", new_callable=PropertyMock
        ) as mock_prop:
            mock_prop.return_value = mock_session

            self.biocyc.clear_caches()
            biocycID = self.biocyc.getDBIdentifierFromInchiKey(
                "GUBGYTABKSRVRQ-ASMJPISFSA-N"
            )
            self.assertEqual("META:ALPHA-MALTOSE", biocycID)

        # Does not exist
        self.biocyc.clear_caches()
        mock_response.text = "Successful queries\nQuery\tBioCyc Common-Name\tBioCyc\tKegg\tCHEBI\tPubChem\tHMDB\tChemSpider\tMetabolights\tMetanetX\tBiGG\tSEED\tInChi\t\n\nAmbigous Queries\nNone\n\nUnknown Queries\nInChIKey:GUBGYTA-KSR-RQ-ASMJ-ISFSA-N\n"
        mock_session.post.return_value = mock_response

        with patch.object(
                Settings, "_biocycSession", new_callable=PropertyMock
        ) as mock_prop:
            mock_prop.return_value = mock_session

            biocycID = self.biocyc.getDBIdentifierFromInchiKey(
                "GUBGYTA-KSR-RQ-ASMJ-ISFSA-N"
            )
            self.assertEqual(Unavailable, biocycID)

        # Returns multiple results
        self.biocyc.clear_caches()
        mock_response.text = 'Successful queries\nQuery\tBioCyc Common-Name\tBioCyc\tKegg\tCHEBI\tPubChem\tHMDB\tChemSpider\tMetabolights\tMetanetX\tBiGG\tSEED\tInChi\t\nInChIKey:GUBGYTABKSRVRQ-DCSYEGIMSA-N\tbeta-lactose\t<A NAME="Sbeta-lactose" HREF="/compound?orgid=META&id=CPD-15971">&beta;-lactose</A>\tC01970\t36218\t6134\tHMDB41627\tNIL\tMTBLC36218\tNIL\tNIL\tNIL\nInChIKey:GUBGYTABKSRVRQ-DCSYEGIMSA-N\tbeta-lactose\t<A NAME="Sbeta-lactose" HREF="/compound?orgid=META&id=CPD-15971">&beta;-lactose</A>\tC01970\t36218\t6134\tHMDB41627\tNIL\tMTBLC36218\tNIL\tNIL\tNIL\n\nAmbigous Queries\nNone\n\nUnknown Queries\nNone\n'
        mock_session.post.return_value = mock_response

        with patch.object(
                Settings, "_biocycSession", new_callable=PropertyMock
        ) as mock_prop:
            mock_prop.return_value = mock_session

            biocycID = self.biocyc.getDBIdentifierFromInchiKey(
                "GUBGYTABKSRVRQ-ASMJPISFSA-N"
            )
            self.assertEqual(Unavailable, biocycID)

        # Live test
        self.biocyc.clear_caches()
        biocycID = self.biocyc.getDBIdentifierFromInchiKey(
            "RGJOEKWQDUBAIZ-IBOSZNHHSA-J"
        )
        self.assertEqual("META:CO-A", biocycID)

        self.biocyc.clear_caches()
        kwargs = {"BioCycSubDB": "ECOLI"}
        biocycID = self.biocyc.getDBIdentifierFromInchiKey(
            "RGJOEKWQDUBAIZ-IBOSZNHHSA-J", **kwargs
        )
        self.assertEqual("ECOLI:CO-A", biocycID)
