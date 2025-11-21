from unittest import TestCase
from unittest.mock import patch, Mock
import requests
from requests import Session

from cobramod.dbwalker.dataclasses import GenerellIdentifiers

from cobramod.dbwalker.BioCyc import BioCyc

import logging

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
        self.biocyc = BioCyc()

    def test_successful_request_with_all_identifiers(self):
        """Test successful API call with all chemical identifiers present."""
        mock_xml = """<ptools-xml ptools-version="29.0" xml:base="http://BioCyc.org/getxml?META:CPD-123">
            <Compound ID="META:CPD-123" orgid="META" frameid="CPD-123" detail="full">
                <cml>
                    <molecule id="CPD-123" title="2,6-dihydroxypyridine" formalCharge="0">
                        <atomArray>
                        </atomArray>
                        <bondArray>
                        </bondArray>
                        <formula concise="C 5 H 5 N 1 O 2"/>
                        <float title="molecularWeight" units="g/mol">111.1</float>
                        <string title="smiles">C1(\C=C(/N=C(\C=1)/O)\O)</string>
                    </molecule>
                </cml>
                <monoisotopic-mw datatype="float" units="Daltons">111.03203</monoisotopic-mw>
                <inchi datatype="string">
                    InChI=1S/C5H5NO2/c7-4-2-1-3-5(8)6-4/h1-3H,(H2,6,7,8)
                </inchi>
                <inchi-key datatype="string">InChIKey=WLFXSECCHULRRO-UHFFFAOYSA-N</inchi-key>

            </Compound>
        </ptools-xml>"""

        with patch.object(Session, "get") as mock_get:
            mock_response = Mock()
            mock_response.content = mock_xml.encode()
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result =  self.biocyc.getGenerellIdentifier("CPD-123")

            self.assertEqual(
                result.inchi, "1S/C5H5NO2/c7-4-2-1-3-5(8)6-4/h1-3H,(H2,6,7,8)"
            )
            self.assertEqual(result.inchi_key, "WLFXSECCHULRRO-UHFFFAOYSA-N")
            self.assertEqual(result.smiles, "C1(\C=C(/N=C(\C=1)/O)\O)")

    def test_successful_request_partial_identifiers(self):
        pass

    def test_custom_database_parameter(self):
        """Test function with custom database parameter."""
        mock_xml = """<?xml version="1.0"?><root></root>"""

        with patch.object(Session, "get") as mock_get:
            mock_response = Mock()
            mock_response.content = mock_xml.encode()
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            self.biocyc.getGenerellIdentifier("CPD-789", BioCycSubDB="ECOLI")

            mock_get.assert_called_once()
            args, kwargs = mock_get.call_args
            assert kwargs["params"]["id"] == "ECOLI:CPD-789"

    def test_request_exception_handling(self):
        """Test handling of requests exceptions."""
        with patch.object(Session, "get") as mock_get:
            mock_get.side_effect = requests.RequestException("Connection error")

            with self.assertLogs(
                "cobramod.DBWalker.BioCyc", level="ERROR"
            ) as cm:
                result =  self.biocyc.getGenerellIdentifier("CPD-123")
                expected = GenerellIdentifiers()

                self.assertEqual(result, expected)
                self.assertEqual(
                    cm.output,
                    [
                        "ERROR:cobramod.DBWalker.BioCyc:Error fetching data from BioCyc: Connection error"
                    ],
                )

    def test_xml_parse_error_handling(self):
        """Test handling of XML parsing errors."""
        invalid_xml = "invalid xml content"

        with patch.object(Session, "get") as mock_get:
            mock_response = Mock()
            mock_response.content = invalid_xml.encode()
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            with self.assertLogs(
                "cobramod.DBWalker.BioCyc", level="ERROR"
            ) as cm:
                result =  self.biocyc.getGenerellIdentifier("CPD-123")
                expected = GenerellIdentifiers()

                self.assertEqual(result, expected)
                expected_error = [
                    "ERROR:cobramod.DBWalker.BioCyc:Error parsing XML response: syntax error: line 1, column 0"
                ]
                self.assertEqual(cm.output, expected_error)

    def test_empty_xml_response(self):
        """Test handling of empty Error site (HTML page)."""
        mock_xml = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
        <html><body>
        </body></html>
        """

        with patch.object(Session, "get") as mock_get:
            mock_response = Mock()
            mock_response.content = mock_xml.encode()
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result =  self.biocyc.getGenerellIdentifier("CPD-999")

            self.assertIsInstance(result, GenerellIdentifiers)
            self.assertEqual(result.inchi, None)
            self.assertEqual(result.inchi_key, None)
            self.assertEqual(result.smiles, None)

    def test_api_parameters(self):
        """Test that correct API parameters are sent."""
        mock_xml = """<?xml version="1.0"?><root></root>"""

        with patch.object(Session, "get") as mock_get:
            mock_response = Mock()
            mock_response.content = mock_xml.encode()
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            self.biocyc.getGenerellIdentifier("CPD-123")

            expected_params = {"id": "META:CPD-123"}
            expected_header = {
                "User-Agent": "Mozilla/5.0",
                "Accept": "application/xml",
            }

            mock_get.assert_called_once_with(
                "https://websvc.biocyc.org/getxml",
                params=expected_params,
                headers=expected_header,
                timeout=30,
            )

    def test_live_api_call(self):
        """Test against the live BioCyc API (may fail if API is unavailable)."""
        result =  self.biocyc.getGenerellIdentifier("GLC")

        self.assertIsInstance(result, GenerellIdentifiers)
        self.assertEqual(
            "C(O)[C@H]1(O[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1)", result.smiles
        )
        self.assertEqual(
            "1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6-/m1/s1",
            result.inchi,
        )
        self.assertEqual(result.inchi_key, "WQZGKKKJIJFFOK-VFUOTHLCSA-N")
